use crate::gf::GF;
use std::ops::{Index, Mul};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

// Macro for creating polynomials, used in tests
#[cfg(test)]
macro_rules! poly {
    ( $( $x:expr ),* ) => {
        {
            let mut temp_vec = Vec::new();
            $(
                temp_vec.push($crate::gf::GF::new($x));
            )*
            $crate::poly::Polynomial::from_slice(&temp_vec)
        }
    };
}

/// Polynomial over GF(2^M) with capacity N coefficients.
///
/// Coefficients are stored in an array `coeffs` where `coeffs[i]`
/// corresponds to the coefficient of x^i. The polynomial degree can
/// be at most N-1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Polynomial<const M: u8, const N: usize> {
    pub coeffs: [GF<M>; N],
}

impl<const M: u8, const N: usize> Polynomial<M, N> {
    /// Creates a polynomial from an array of coefficients.
    pub fn new(coeffs: [GF<M>; N]) -> Self {
        Polynomial { coeffs }
    }

    /// Returns the zero polynomial (all coefficients set to 0).
    pub const fn zero() -> Self {
        Polynomial { coeffs: [GF(0); N] }
    }

    /// Creates a polynomial from a slice of field elements.
    ///
    /// If the slice is shorter than N, the remaining coefficients
    /// are set to zero. If the slice is longer than N, it is truncated.
    pub fn from_slice(slice: &[GF<M>]) -> Self {
        let mut coeffs = [GF::new(0); N];

        let copy_len = slice.len().min(N);
        coeffs[..copy_len].copy_from_slice(&slice[..copy_len]);

        Polynomial { coeffs }
    }

    /// Computes the degree of the polynomial in constant time.
    ///
    /// Returns the largest index i such that coeffs[i] != 0, or 0
    /// if the polynomial is the zero polynomial.
    ///
    /// # Constant-time
    /// Yes. Scans all N coefficients using bitwise operations and
    /// conditional selection.
    pub fn deg(&self) -> usize {
        let mut res = 0u32;
        let mut found = Choice::from(0);
        for i in (0..N).rev() {
            let is_nonzero = !self.coeffs[i].ct_eq(&GF::new(0));
            let update = is_nonzero & !found;
            res = u32::conditional_select(&res, &(i as u32), update);
            found |= is_nonzero;
        }
        res as usize
    }

    /// Returns the leading coefficient of the polynomial.
    ///
    /// The leading coefficient is the coefficient of the highest-degree
    /// non-zero term. Returns 0 if the polynomial is zero.
    ///
    /// # Constant-time
    /// Yes. Scans all N coefficients using conditional selection.
    pub fn lead(&self) -> GF<M> {
        let mut leading = GF::new(0);
        let mut found = Choice::from(0);
        for i in (0..N).rev() {
            let is_nonzero = !self.coeffs[i].ct_eq(&GF::new(0));
            let update = is_nonzero & !found;

            leading = GF::conditional_select(&leading, &self.coeffs[i], update);
            found |= is_nonzero;
        }
        leading
    }

    /// Checks whether the polynomial is the zero polynomial.
    ///
    /// Returns `Choice(1)` if all coefficients are zero, `Choice(0)`
    /// otherwise.
    ///
    /// # Constant-time
    /// Yes. ANDs together the equality checks for all coefficients.
    pub fn is_zero(&self) -> Choice {
        let mut res = Choice::from(1);
        for i in 0..N {
            res &= self.coeffs[i].ct_eq(&GF::new(0));
        }
        res
    }

    /// Conditionally swaps two polynomials in constant time.
    ///
    /// If `choice` is `Choice(1)`, `a` and `b` are swapped; otherwise
    /// they remain unchanged.
    ///
    /// # Constant-time
    /// Yes. Performs a masked swap on every coefficient pair.
    #[inline(always)]
    pub fn swap(a: &mut Self, b: &mut Self, choice: Choice) {
        for i in 0..N {
            let temp_a = GF::conditional_select(&a.coeffs[i], &b.coeffs[i], choice);
            let temp_b = GF::conditional_select(&b.coeffs[i], &a.coeffs[i], choice);
            a.coeffs[i] = temp_a;
            b.coeffs[i] = temp_b;
        }
    }

    /// Shifts all coefficients down by one position.
    ///
    /// This corresponds to dividing the polynomial by x.
    /// The last (highest-index) coefficient is set to zero.
    ///
    /// # Constant-time
    /// Yes. Iterates over all N-1 coefficients unconditionally.
    #[inline(always)]
    pub fn shift_down(&mut self) {
        for i in 0..N - 1 {
            self.coeffs[i] = self.coeffs[i + 1];
        }
        self.coeffs[N - 1] = GF::new(0);
    }

    /// Evaluates the polynomial at `x` using Horner's method.
    ///
    /// # Constant-time
    /// Yes. Iterates over all N coefficients unconditionally.
    pub fn eval(&self, x: GF<M>) -> GF<M> {
        let mut res = GF::new(0);
        for i in (0..N).rev() {
            res = (res * x) + self.coeffs[i];
        }
        res
    }

    /// Normalizes the polynomial to monic form in constant time.
    ///
    /// Multiplies every coefficient by the inverse of the leading
    /// coefficient. If the polynomial is the zero polynomial, the
    /// operation is masked to avoid division by zero.
    ///
    /// # Constant-time
    /// Yes. Uses conditional selection throughout to avoid leaking
    /// the position or value of the leading coefficient.
    pub fn make_monic(&mut self) {
        let mut leading = GF::new(0);
        let mut found = Choice::from(0);

        // Find the leading coefficient (without using array indexing that cause cache-timing attacks)
        for i in (0..N).rev() {
            let is_nonzero = !self.coeffs[i].ct_eq(&GF::new(0));
            let update = is_nonzero & !found;
            leading = GF::conditional_select(&leading, &self.coeffs[i], update);
            found |= is_nonzero;
        }

        // If the leading coefficient is zero, mask it to avoid division by zero
        let is_poly_zero = !found;
        let safe_leading = GF::conditional_select(&leading, &GF::new(1), is_poly_zero);
        let inv = safe_leading.inv();

        // Multiply each coefficient by the inverse of the leading coefficient to normalize the polynomial
        for i in 0..N {
            self.coeffs[i] = self.coeffs[i] * inv;
        }
    }

    /// Constant-time polynomial division with remainder.
    ///
    /// Computes `dividend / divisor` and `dividend % divisor` where
    /// the divisor must be monic and non-zero. `d` is the degree of
    /// the divisor.
    ///
    /// # Constant-time
    /// Yes. All loops iterate a fixed number of times (N - d iterations
    /// for the outer loop, d iterations for the inner loop).
    pub fn ct_div_rem(dividend: &Self, divisor: &Self, d: usize) -> (Self, Self) {
        debug_assert!(
            divisor.deg() == d,
            "ct_div_rem: d ({}) must equal divisor degree ({})",
            d,
            divisor.deg()
        );
        debug_assert!(
            divisor.is_zero().unwrap_u8() == 0,
            "ct_div_rem: divisor must be non-zero"
        );

        let mut rem = *dividend;
        let mut q = Self::zero();

        for i in (d..N).rev() {
            let coef = rem.coeffs[i];
            q.coeffs[i - d] = coef;

            // Subtract coef * divisor * x^(i-d) from the remainder
            for j in 0..d {
                rem.coeffs[i - d + j] = rem.coeffs[i - d + j] + (divisor.coeffs[j] * coef);
            }

            rem.coeffs[i] = GF::new(0);
        }

        (q, rem)
    }

    /// Bernstein-Yang SafeGCD algorithm for constant-time polynomial GCD.
    ///
    /// Computes `gcd(a, b)` over GF(2^M)[x]. The result is normalized
    /// to monic form.
    ///
    /// # Constant-time
    /// Yes. Runs for exactly `2 * max_deg + 1` iterations regardless of
    /// the input polynomials. All operations use conditional selection.
    pub fn gcd(a: &Self, b: &Self, max_deg: usize) -> Self {
        let mut f = *a;
        let mut g = *b;
        let mut delta = 1i32;

        for _ in 0..=(2 * max_deg) {
            let delta_gt_0 = Choice::from(u8::from(delta > 0));
            let g0_not_0 = !g.coeffs[0].ct_eq(&GF::new(0));
            let cond = delta_gt_0 & g0_not_0;

            let neg_delta = -delta;
            delta = i32::conditional_select(&delta, &neg_delta, cond) + 1;

            Self::swap(&mut f, &mut g, cond);

            let f0 = f.coeffs[0];
            let g0 = g.coeffs[0];

            let mut g_new = Self::zero();
            for i in 0..N {
                g_new.coeffs[i] = (f0 * g.coeffs[i]) + (g0 * f.coeffs[i]);
            }
            g_new.shift_down();

            g = g_new;
        }
        f.make_monic();

        f
    }

    /// Repeated square-and-multiply for polynomial exponentiation in F_q^m.
    ///
    /// Computes `self^k mod f` where `f` has degree `f_deg`.
    ///
    /// The iteration count is fixed at the bit-width of `usize` (typically 64),
    /// making this constant-time regardless of the exponent's Hamming weight.
    /// The result is correct only when k < 2^64 — which holds for all internal
    /// callers (k = 1 << M with M ∈ {12, 13}).
    ///
    /// Reference: Handbook of Applied Cryptography, Algorithm 2.227
    ///
    /// # Constant-time
    /// Yes. The loop iterates exactly `usize::BITS` times regardless of
    /// the exponent value.
    pub fn mod_pow(&self, k: usize, f: &Self, f_deg: usize) -> Self {
        let mut s = Self::zero();
        s.coeffs[0] = GF(1);

        let mut g_x = self.reduce(f, f_deg);

        let mut current_k = k;

        for _ in 0..(usize::BITS as usize) {
            let bit = Choice::from((current_k & 1) as u8);

            let prod = &s * &g_x;
            let reduced = prod.reduce(f, f_deg);
            for i in 0..N {
                s.coeffs[i] = GF::conditional_select(&s.coeffs[i], &reduced.coeffs[i], bit);
            }

            let sq = &g_x * &g_x;
            g_x = sq.reduce(f, f_deg);

            current_k >>= 1;
        }

        s
    }

    /// Ben-Or irreducibility test (constant-time variant).
    ///
    /// Tests whether the polynomial is irreducible over GF(2^M)[x].
    ///
    /// Reference: Handbook of Applied Cryptography, Algorithm 4.69
    ///
    /// # Constant-time
    /// Yes. The algorithm has been modified from the original to use
    /// fixed iteration counts and avoid early termination.
    pub fn is_irreducible(&self, expected_deg: usize) -> Choice {
        let mut u = Self::zero();
        u.coeffs[1] = GF::new(1);
        let q = 1usize << (M as usize);

        let mut is_irred = Choice::from(1);

        for _ in 1..=(expected_deg / 2) {
            u = u.mod_pow(q, self, expected_deg);

            let mut u_minus_x = u;
            u_minus_x.coeffs[1] = u_minus_x.coeffs[1] + GF::new(1);

            let d = Self::gcd(&u_minus_x, self, expected_deg);

            let deg_is_zero = Choice::from((d.deg() == 0) as u8);
            is_irred &= deg_is_zero;
        }

        let correct_deg = Choice::from((self.deg() == expected_deg) as u8);

        is_irred & correct_deg
    }

    /// Reduction of a polynomial modulo `f`.
    ///
    /// Computes `self mod f` where `f` has degree `t`.
    ///
    /// # Constant-time
    /// Yes. Iterates from degree t to N-1 unconditionally with
    /// branchless operations.
    pub fn reduce(&self, f: &Self, t: usize) -> Self {
        let mut r = self.coeffs;
        for i in (t..N).rev() {
            let c = r[i];
            for j in 0..t {
                let f_j = f.coeffs[j];
                r[i - t + j] = r[i - t + j] + (c * f_j);
            }
            r[i] = GF::new(0);
        }
        Polynomial::new(r)
    }

    /// Divide-and-conquer product tree for polynomial multiplication.
    ///
    /// Recursively multiplies a list of polynomials together,
    /// reducing modulo `f_y` at each step.
    #[allow(dead_code)]
    fn product_tree(factors: &[Self], f_y: &Self) -> Self {
        match factors.len() {
            0 => {
                let mut p = Self::zero();
                p.coeffs[0] = GF::new(1);
                p
            } // empty product = 1
            1 => factors[0], // base case
            _ => {
                // split down the middle
                let mid = factors.len() / 2;
                let left = Self::product_tree(&factors[..mid], f_y);
                let right = Self::product_tree(&factors[mid..], f_y);
                // multiply and reduce mod f_y
                let prod = &left * &right;
                prod.reduce(f_y, f_y.deg())
            }
        }
    }

    /// Extended product tree building outer polynomial coefficients.
    ///
    /// Given conjugates [a_0, ..., a_{k-1}], builds the outer polynomial
    /// (X - a_0)(X - a_1)...(X - a_{k-1}) via divide-and-conquer.
    fn product_tree_ext(factors: &[Self], f_y: &Self) -> Vec<Self> {
        match factors.len() {
            0 => {
                let mut p = Self::zero();
                p.coeffs[0] = GF::new(1);
                vec![p]
            }
            1 => {
                // (X + conj) = [conj, 1]
                let mut p = Self::zero();
                p.coeffs[0] = GF::new(1);
                vec![factors[0], p]
            }
            _ => {
                let mid = factors.len() / 2;
                let left = Self::product_tree_ext(&factors[..mid], f_y);
                let right = Self::product_tree_ext(&factors[mid..], f_y);

                // multiply left and right as outer polynomials
                // reuse your existing loop logic, just extracted here
                let mut res = vec![Self::zero(); left.len() + right.len() - 1];
                for (i, ca) in left.iter().enumerate() {
                    for (j, cb) in right.iter().enumerate() {
                        let prod = ca * cb;
                        let rem = Self::reduce(&prod, f_y, f_y.deg());

                        let r = &mut res[i + j];

                        for k in 0..N {
                            r.coeffs[k] = r.coeffs[k] + rem.coeffs[k];
                        }
                    }
                }
                res
            }
        }
    }

    /// Frobenius squaring in characteristic 2.
    ///
    /// Computes `p^2 mod f`. In characteristic 2, (sum a_i y^i)^2 =
    /// sum a_i^2 y^(2i) since cross terms vanish.
    fn frobenius_sq(p: &Self, f: &Self) -> Self {
        let deg = p.deg();
        let mut res = Self::zero();
        for i in 0..=deg {
            res.coeffs[2 * i] = p.coeffs[i].sq(); // squaring each GF element, coefficients go to even positions
        }
        res.reduce(f, f.deg())
    }

    /// Full Frobenius automorphism: p -> p^(2^M) mod f.
    ///
    /// Instead of `mod_pow(2^M, f)` which does 2^M iterations,
    /// we do M squarings — from 8192 iterations down to 13.
    fn frobenius(&self, f: &Self) -> Self {
        let mut result = *self;
        for _ in 0..M {
            result = Polynomial::frobenius_sq(&result, f);
        }
        result
    }

    /// Computes the minimal polynomial of `self` in GF(2^M)[y].
    ///
    /// Collects conjugates via the Frobenius automorphism and
    /// multiplies them as (X - conj) factors.
    ///
    /// # Constant-time
    /// Partially. The number of conjugates is fixed at `expected_deg`,
    /// but the result uses heap allocation for the product tree.
    pub fn minpoly(&self, f_y: &Self, expected_deg: usize) -> Self {
        let mut conjugates: Vec<Self> = Vec::with_capacity(expected_deg);
        let mut current = *self;

        for _ in 0..expected_deg {
            conjugates.push(current);
            current = current.frobenius(f_y);
        }

        // Multiply out (X - conj_0)(X - conj_1)...
        let acc = Polynomial::product_tree_ext(&conjugates, f_y);

        // At this point acc[i] should each be a degree-0 polynomial (a scalar in GF<M>)
        // because minpoly lands back in GF(2^M)[y] — extract those scalars
        let mut result = Self::zero();
        let len = acc.len().min(N);
        for (i, item) in acc.iter().enumerate().take(len) {
            result.coeffs[i] = item.coeffs[0];
        }

        result.make_monic();
        result
    }

    /// Bernstein-Yang SafeGCD for polynomial inversion modulo `m`.
    ///
    /// Returns `(inverse, is_invertible)` where `is_invertible` is
    /// `Choice(1)` if the inverse exists and `Choice(0)` otherwise.
    ///
    /// # Constant-time
    /// Yes. Runs for exactly `2 * max_deg + 1` iterations. All
    /// data-dependent branches use conditional selection.
    pub fn inv_mod(&self, m: &Self, max_deg: usize) -> (Self, Choice) {
        let mut f = *m;
        let mut g = *self;

        let mut v = Self::zero();
        let mut r = Self::zero();
        r.coeffs[0] = GF::new(1);

        let mut delta = 1i32;

        let m0_inv = m.coeffs[0].inv();

        for _ in 0..=(2 * max_deg) {
            let delta_gt_0 = Choice::from(u8::from(delta > 0));
            let g0_not_zero = !g.coeffs[0].ct_eq(&GF::new(0));
            let cond = delta_gt_0 & g0_not_zero;

            let neg_delta = -delta;
            delta = i32::conditional_select(&delta, &neg_delta, cond) + 1;

            Self::swap(&mut f, &mut g, cond);
            Self::swap(&mut v, &mut r, cond);

            let f0 = f.coeffs[0];
            let g0 = g.coeffs[0];

            let mut g_new = Self::zero();
            let mut r_new = Self::zero();
            for i in 0..N {
                g_new.coeffs[i] = (f0 * g.coeffs[i]) + (g0 * f.coeffs[i]);
                r_new.coeffs[i] = (f0 * r.coeffs[i]) + (g0 * v.coeffs[i]);
            }
            g_new.shift_down();

            let c = r_new.coeffs[0] * m0_inv;
            for i in 0..N {
                r_new.coeffs[i] = (c * m.coeffs[i]) + r_new.coeffs[i];
            }
            r_new.shift_down();

            g = g_new;
            r = r_new;
        }

        let f0_final = f.coeffs[0];
        let is_invertible = !f0_final.ct_eq(&GF::new(0));

        let safe_f0 = GF::conditional_select(&f0_final, &GF::new(1), !is_invertible);
        let safe_f0_inv = safe_f0.inv();

        let mut inverse = Self::zero();
        for i in 0..N {
            inverse.coeffs[i] = v.coeffs[i] * safe_f0_inv;
        }

        (inverse, is_invertible)
    }
}

/// Index into a polynomial's coefficient array.
impl<const M: u8, const N: usize> Index<usize> for Polynomial<M, N> {
    type Output = GF<M>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coeffs[index]
    }
}

/// Polynomial multiplication (convolution) in GF(2^M)[x].
///
/// Multiplies two polynomials using schoolbook O(N^2) multiplication.
/// The result is truncated to capacity N (coefficients of degree >= N
/// are silently dropped).
///
/// # Constant-time
/// Yes. All loops iterate N×N times unconditionally.
impl<const M: u8, const N: usize> Mul for &Polynomial<M, N> {
    type Output = Polynomial<M, N>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = Polynomial::<M, N>::zero();

        for i in 0..N {
            for j in 0..N {
                if i + j < N {
                    res.coeffs[i + j] = res.coeffs[i + j] + (self[i] * rhs[j]);
                }
            }
        }

        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf::GF;

    type TestGF = GF<13>;
    type TestPoly = Polynomial<13, 256>;
    const T: usize = 96;

    fn f_y_460896() -> TestPoly {
        // F(y) = y^96 + y^10 + y^9 + y^6 + 1, z = GF<13>(2)
        let mut test_poly = TestPoly::zero();
        test_poly.coeffs[0] = TestGF::new(1); // 1
        test_poly.coeffs[6] = TestGF::new(1); // y^6
        test_poly.coeffs[9] = TestGF::new(1); // y^9
        test_poly.coeffs[10] = TestGF::new(1); // y^10
        test_poly.coeffs[96] = TestGF::new(1); // y^96
        test_poly
    }

    #[test]
    fn test_poly_eval() {
        let p: TestPoly = poly![1, 1, 1]; // x^2 + x + 1
        let x_val = TestGF::new(2);
        let result = p.eval(x_val);
        assert_eq!(result, TestGF::new(7));
    }

    #[test]
    fn test_poly_ct_div_rem() {
        let p: TestPoly = poly![1, 1, 1, 1]; // x^3 + x^2 + x + 1
        let q: TestPoly = poly![1, 1]; // x + 1
        let (div, rem) = TestPoly::ct_div_rem(&p, &q, q.deg());
        assert_eq!(div, poly![1, 0, 1]); // x^2 + 1
        assert_eq!(rem.is_zero().unwrap_u8(), 1);
    }

    #[test]
    fn test_poly_gcd() {
        let p: TestPoly = poly![1, 0, 0, 1]; // x^3 + 1
        let q: TestPoly = poly![1, 0, 1]; // x^2 + 1
        let gcd = TestPoly::gcd(&p, &q, p.deg());
        assert_eq!(gcd, poly![1, 1]); // x + 1
    }

    #[test]
    fn test_poly_mod_pow() {
        let p: TestPoly = poly![0, 1]; // x
        let q: TestPoly = poly![1, 1, 1]; // x^2 + x + 1
        let result = p.mod_pow(3, &q, q.deg());
        assert_eq!(result, poly![1]); //  1
        assert_eq!(result.deg(), 0);
    }

    #[test]
    fn test_poly_is_irreducible() {
        let p: TestPoly = poly![1, 1, 0, 1]; // x^3 + x + 1
        assert!(p.is_irreducible(3).unwrap_u8() == 1);
        let q: TestPoly = poly![1, 0, 0, 1]; // x^3 + 1 = (x + 1)(x^2 + x + 1)
        assert!(q.is_irreducible(3).unwrap_u8() == 0);
    }

    #[test]
    fn test_poly_minpoly() {
        let f_y = f_y_460896();

        // beta = 1 + y + y^2
        let mut beta = TestPoly::zero();
        beta.coeffs[0] = TestGF::new(1);
        beta.coeffs[1] = TestGF::new(1);
        beta.coeffs[2] = TestGF::new(1);

        let g = beta.minpoly(&f_y, f_y.deg());

        assert_eq!(g.coeffs[g.deg()].0, 1, "minpoly must be monic");
        assert!(
            T % g.deg() == 0,
            "deg(g) = {} must divide T = {}",
            g.deg(),
            T
        );
        assert!(
            g.is_irreducible(g.deg()).unwrap_u8() == 1,
            "minpoly must be irreducible"
        );

        // g is USE 1: a real polynomial with GF13 scalar coefficients
        // beta is USE 2: an extension ring element represented as Polynomial<13>

        // g(beta) means: for each term g[i] * y^i in g,
        //   substitute beta for y -> g[i] * beta^i
        //   where g[i] is a GF13 scalar  (scales the extension element)
        //   and   beta^i is computed via repeated mul + div_rem mod f_y
        let mut result = TestPoly::zero();
        let mut beta_pow = TestPoly::zero(); // beta^0 = 1
        beta_pow.coeffs[0] = TestGF::new(1);

        for i in 0..=g.deg() {
            let scaled: Vec<TestGF> = beta_pow.coeffs.iter().map(|&c| c * g.coeffs[i]).collect();
            let scaled_poly = TestPoly::from_slice(&scaled);

            let len = result.coeffs.len().max(scaled_poly.coeffs.len());
            let mut res = vec![TestGF::new(0); len];
            for (j, c) in result.coeffs.iter().enumerate() {
                res[j] = res[j] + *c;
            }
            for (j, c) in scaled_poly.coeffs.iter().enumerate() {
                res[j] = res[j] + *c;
            }
            result = TestPoly::from_slice(&res);

            let prod = &beta_pow * &beta;
            beta_pow = prod.reduce(&f_y, f_y.deg());
        }

        assert_eq!(
            result.is_zero().unwrap_u8(),
            1,
            "g(beta) must be 0 in GF(2^13), got {:?}",
            result
        );
    }

    #[test]
    fn test_poly_inv_mod() {
        let f = f_y_460896();

        // a(x) = x^2 + x + 1
        let a = poly![1, 1, 1];

        let (inv, inv_valid) = a.inv_mod(&f, f.deg());
        assert!(
            inv_valid.unwrap_u8() == 1,
            "a(x) must have an inverse modulo f(x)"
        );

        // a(x) * a^-1(x)
        let prod = &a * &inv;

        let rem = TestPoly::reduce(&prod, &f, f.deg());

        assert_eq!(rem.deg(), 0, "a(x) * a^-1(x) mod f(x) should have degree 0");
        assert_eq!(rem.coeffs[0].0, 1, "a(x) * a^-1(x) mod f(x) must be 1");
    }
}
