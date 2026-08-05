use std::ops::{Add, Mul};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// Element of the finite field GF(2^M).
///
/// Internally represented as a polynomial in GF(2)[x] modulo the
/// irreducible polynomial of degree M, stored in a `u16`.
///
/// Implements `ConditionallySelectable` and `ConstantTimeEq` for
/// use in constant-time cryptographic operations.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct GF<const M: u8>(pub u16);

impl<const M: u8> GF<M> {
    /// Creates a new field element, reducing the value modulo 2^M.
    ///
    /// Values are masked to ensure they fit within the field size.
    pub fn new(value: u16) -> Self {
        GF(value & ((1u16 << M) - 1))
    }

    /// Returns the zero element of the field (all bits set to 0).
    pub fn zero() -> Self {
        GF(0)
    }

    /// Returns the irreducible polynomial for GF(2^M) as a bitmask.
    ///
    /// Used internally for carry-less multiplication reduction.
    /// Currently supports M = 12 (0x1009) and M = 13 (0x201B).
    pub fn get_irreducible_poly() -> u32 {
        match M {
            12 => 0x1009,
            13 => 0x201B,
            _ => panic!("No irreducible defined!"),
        }
    }

    /// Fast inverse using Fermat's little theorem: a^(2^M - 2) mod p.
    ///
    /// Zero maps to zero (0^-1 = 0).
    ///
    /// # Constant-time
    /// Yes. Delegates to [`pow`](Self::pow) which uses a fixed iteration count.
    pub fn inv(self) -> Self {
        self.pow((1 << M) - 2)
    }

    /// Squaring in GF(2^M).
    ///
    /// Computes `self * self`. In characteristic 2, squaring is linear:
    /// cross terms vanish, making the result computable via `mul`.
    ///
    /// # Constant-time
    /// Yes. Delegates to [`mul`](Mul::mul) which is constant-time.
    ///
    /// TODO: Use bit interleaving for faster/optimized squaring
    pub fn sq(self) -> Self {
        self.mul(self)
    }

    /// Exponentiation in GF(2^M) via square-and-multiply.
    ///
    /// Computes `self^exp` where `exp` is a 16-bit exponent.
    ///
    /// # Constant-time
    /// Yes. The loop iterates exactly 16 times regardless of the
    /// exponent value or Hamming weight.
    pub fn pow(self, exp: u16) -> Self {
        let mut res = GF::<M>::new(1);
        let mut base = self;
        let mut current_exp = exp;

        for _ in 0..16 {
            let bit = (current_exp & 1) as u8;
            let choice = Choice::from(bit);
            let prod = res.mul(base);
            res = GF::conditional_select(&res, &prod, choice);
            base = base.sq();
            current_exp >>= 1;
        }

        res
    }

    /// Square root in GF(2^M).
    ///
    /// Computes `a^(2^(M-1))`. In characteristic 2, squaring is a
    /// bijection, so every element has a unique square root.
    ///
    /// # Constant-time
    /// Yes. Delegates to [`pow`](Self::pow) which is constant-time.
    pub fn sqrt(self) -> Self {
        self.pow(1 << (M - 1))
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
/// Addition in GF(2^M) — bitwise XOR.
///
/// # Constant-time
/// Yes. XOR is a basic arithmetic operation with no data-dependent
/// timing.
impl<const M: u8> Add for GF<M> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        GF(self.0 ^ rhs.0)
    }
}

/// Carry-less multiplication in GF(2^M).
///
/// Multiplies two field elements using polynomial multiplication
/// followed by reduction modulo the irreducible polynomial.
///
/// # Constant-time
/// Yes. All loops iterate a fixed number of times (M iterations for
/// the product, M iterations for the reduction) and use branchless
/// masks for conditional operations.
impl<const M: u8> Mul for GF<M> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut res: u32 = 0u32;
        let m_mask = (1 << M) - 1;

        // Cast as u32 to avoid overflow
        // & m_mask to ensure the value is within the range [0, 2^M - 1]
        let a = (self.0 as u32) & m_mask;
        let b = (rhs.0 as u32) & m_mask;

        // Polynomial multiplication a(x)b(x)
        for i in 0..M {
            let bit = (a >> i) & 1;
            // let mask = -(bit as i16) as u32;
            // Branchless mask: 0xfffffffe if bit == 0, 0x00000000 if bit == 1
            let mask = 0u32.wrapping_sub(bit);
            res ^= (b << i) & mask;
        }

        // Polynomial reduction mod p(x)
        let poly: u32 = Self::get_irreducible_poly();
        for i in (M..=(2 * M - 1)).rev() {
            let bit = (res >> i) & 1;
            let mask = 0u32.wrapping_sub(bit);
            res ^= (poly << (i - M)) & mask;
        }

        GF((res & m_mask) as u16)
    }
}

/// Constant-time equality for GF elements.
///
/// # Constant-time
/// Yes. Delegates to the underlying `u16` constant-time equality.
impl<const M: u8> ConstantTimeEq for GF<M> {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

/// Constant-time conditional selection for GF elements.
///
/// # Constant-time
/// Yes. Delegates to `u16::conditional_select`.
impl<const M: u8> ConditionallySelectable for GF<M> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        GF(u16::conditional_select(&a.0, &b.0, choice))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    type TestGF = GF<13>;

    #[test]
    fn test_gf_new() {
        // Test that values are reduced modulo the irreducible polynomial
        let a = TestGF::new(8192);
        assert_eq!(a.0, 0);

        let b = TestGF::new(8195);
        assert_eq!(b.0, 3);
    }

    #[test]
    fn test_gf_add() {
        let a = TestGF::new(5);
        let b = TestGF::new(3);
        let zero = TestGF::new(0);
        assert_eq!(a + a, zero);
        assert_eq!(a + zero, a);
        assert_eq!(a + b, TestGF::new(6));
        assert_eq!(a + b, b + a);
    }

    #[test]
    fn test_gf_mul() {
        let a = TestGF::new(5);
        let b = TestGF::new(3);
        let zero = TestGF::new(0);
        let one = TestGF::new(1);
        assert_eq!(a * a, TestGF::new(17)); // (x^2 + 1)(x^2 + 1) = x^4 + 1
        assert_eq!(a * zero, zero);
        assert_eq!(a * b, TestGF::new(15));
        assert_eq!(a * b, b * a);
        assert_eq!(a * a * a, TestGF::new(85));
        assert_eq!(a * one, a);
        assert_eq!(one * a, a);
    }

    #[test]
    fn test_gf_sq() {
        let a = TestGF::new(5);
        assert_eq!(a.sq(), TestGF::new(17)); // (x^2 + 1)^2 = x^4 + 1
    }

    #[test]
    fn test_gf_pow() {
        let a = TestGF::new(2);
        assert_eq!(a.pow(3), TestGF::new(8));
        assert_eq!(a.pow(0), TestGF::new(1));
    }

    #[test]
    fn test_gf_inv() {
        let a = TestGF::new(15);
        let a_inv = a.inv();
        assert_eq!(a_inv * a, TestGF::new(1)); // 15 * 11 = 1
    }

    #[test]
    fn test_gf_sqrt() {
        let a = TestGF::new(1234);
        let a_sq = a.sq();

        assert_eq!(a_sq.sqrt(), a, "Sqrt failed");

        let zero = TestGF::new(0);
        assert_eq!(zero.sqrt(), zero, "sqrt(0) should be 0");

        let one = TestGF::new(1);
        assert_eq!(one.sqrt(), one, "sqrt(1) should be 1");
    }

    /// Uniform strategy over all 8192 elements of GF(2^13): 0..2^M.
    ///
    /// A local strategy beats `impl Arbitrary for GF<M>` here: GF is
    /// generic over `const M`, so an `Arbitrary` impl would have to be
    /// a blanket impl for all M, leaking proptest into the public API
    /// (or requiring `#[cfg(test)]` gating on a trait impl).
    fn any_gf13() -> impl Strategy<Value = TestGF> {
        (0u16..(1 << 13)).prop_map(TestGF::new)
    }

    proptest! {
        /// Addition in GF(2^m) is XOR, so (GF, +) is an abelian group:
        /// (a + b) + c == a + (b + c) must hold for every a, b, c.
        #[test]
        fn prop_gf_add_associative(a in any_gf13(), b in any_gf13(), c in any_gf13()) {
            prop_assert_eq!((a + b) + c, a + (b + c));
        }

        /// Multiplication in GF(2^m) is defined by the polynomial modulus, so (GF, *) is a group:
        /// (a * b) * c == a * (b * c) must hold for every a, b, c.
        #[test]
        fn prop_gf_mul_associative(a in any_gf13(), b in any_gf13(), c in any_gf13()) {
            prop_assert_eq!((a * b) * c, a * (b * c));
        }

        /// Addition in GF(2^m) is XOR, so (GF, +) is an abelian group:
        /// a + b == b + a must hold for every a, b.
        #[test]
        fn prop_gf_add_commutative(a in any_gf13(), b in any_gf13()) {
            prop_assert_eq!(a + b, b + a);
        }

        /// Multiplication in GF(2^m) is defined by the polynomial modulus, so (GF, *) is a group:
        /// a * b == b * a must hold for every a, b.
        #[test]
        fn prop_gf_mul_commutative(a in any_gf13(), b in any_gf13()) {
            prop_assert_eq!(a * b, b * a);
        }

        /// Multiplication in GF(2^m) is distributive:
        /// a * (b + c) == a * b + a * c must hold for every a, b, c.
        #[test]
        fn prop_gf_mul_distributive(a in any_gf13(), b in any_gf13(), c in any_gf13()) {
            prop_assert_eq!(a * (b + c), a * b + a * c);
        }
    }
}
