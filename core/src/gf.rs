use std::ops::{Add, Mul};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct GF<const M: u8>(pub u16);

impl<const M: u8> GF<M> {
    pub fn new(value: u16) -> Self {
        GF(value & ((1u16 << M) - 1))
    }

    pub fn get_irreducible_poly() -> u32 {
        match M {
            12 => 0x1009,
            13 => 0x201B,
            _ => panic!("No irreducible defined!"),
        }
    }

    // Fast inverse using fermat's little theorem (a^(2^M-2) mod p) (0^-1 is 0)
    pub fn inv(self) -> Self {
        self.pow((1 << M) - 2)
    }

    // TODO: Use bit interleaving for faster/optimized squaring
    pub fn sq(self) -> Self {
        self.mul(self)
    }

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

    pub fn sqrt(self) -> Self {
        self.pow(1 << (M - 1))
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const M: u8> Add for GF<M> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        GF(self.0 ^ rhs.0)
    }
}

// Carry-less multiplication
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

impl<const M: u8> ConstantTimeEq for GF<M> {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl<const M: u8> ConditionallySelectable for GF<M> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        GF(u16::conditional_select(&a.0, &b.0, choice))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
