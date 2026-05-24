use crate::gf::GF;
use std::ops::{Index, Mul};

// Macro for creating polynomials
macro_rules! poly {
    ( $( $x:expr ),* ) => {
        {
            let mut temp_vec = Vec::new();
            $(
                temp_vec.push(crate::gf::GF::new($x));
            )*
            crate::poly::Polynomial::new(temp_vec)
        }
    };
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial<const M: u8> {
    pub coeffs: Vec<GF<M>>,
}

impl<const M: u8> Polynomial<M> {
    pub fn new(coeffs: Vec<GF<M>>) -> Self {
        Polynomial { coeffs }
    }

    pub fn clean(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last().unwrap().0 == 0 {
            self.coeffs.pop();
        }
    }

    pub fn deg(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || (self.coeffs.len() == 1 && self.coeffs[0].0 == 0)
    }

    // Horner's method for polynomial evaluation
    pub fn eval(&self, x: GF<M>) -> GF<M> {
        let mut res = GF::new(0);
        for &coeff in self.coeffs.iter().rev() {
            res = (res * x) + coeff;
        }
        res
    }

    pub fn make_monic(&mut self) {
        if self.is_zero() {
            return;
        }
        let leading = *self.coeffs.last().unwrap();
        if leading.0 != 1 {
            let inv = leading.inv();
            for coeff in self.coeffs.iter_mut() {
                *coeff = *coeff * inv;
            }
        }
    }

    pub fn div_rem(dividend: &Self, divisor: &Self) -> (Self, Self) {
        let mut rem = dividend.clone();
        rem.clean();

        let mut div = divisor.clone();
        div.clean();

        if div.is_zero() {
            panic!("divisor is zero");
        }

        // If the degree of the remainder is less than the degree of the divisor, return 0 and the remainder
        if rem.deg() < div.deg() {
            return (Polynomial::new(vec![GF::new(0)]), rem);
        }

        let mut q_coeffs = vec![GF::new(0); dividend.deg() - divisor.deg() + 1];

        let div_lead_inv = div.coeffs.last().unwrap().inv();

        while !rem.is_zero() && rem.deg() >= div.deg() {
            let deg_diff = rem.deg() - div.deg();
            let rem_lead = *rem.coeffs.last().unwrap();

            let ratio = rem_lead * div_lead_inv;

            q_coeffs[deg_diff] = ratio;

            // Remainder = Remainder - (ration * divisor) [in XOR (+) and (-) is the same]
            for i in 0..div.coeffs.len() {
                rem.coeffs[deg_diff + i] = rem.coeffs[deg_diff + i] + (div.coeffs[i] * ratio);
            }

            rem.clean();
        }

        let mut q = Polynomial::new(q_coeffs);
        q.clean();

        (q, rem)
    }

    // Handbook of Applied Cryptography, Algorithm 2.218, Euclidean Algorithm for Z_p[x]
    // INPUT: two polynomials g and h over Z_p[x]
    // OUTPUT: the greatest common divisor of g and h
    pub fn gcd(g: &Self, h: &Self) -> Self {
        let mut g = g.clone();
        let mut h = h.clone();
        while !h.is_zero() {
            let (_, r) = Self::div_rem(&g, &h);
            g = h;
            h = r;
        }
        g.make_monic();
        g
    }

    // Handbook of Applied Cryptography, Algorithm 2.227, Repeated square-and-multiply algorithm for exponentiation in F_q^m
    // INPUT: a polynomial g in F_q^m (&self), and an integer 0 <= k < p^m - 1 (where F_q^m = Z_p[x]/f)
    // OUTPUT: the result of g^k mod f
    pub fn mod_pow(&self, mut k: usize, f: &Self) -> Self {
        let mut s = Polynomial::new(vec![GF::new(1)]);
        if k == 0 {
            return s;
        }

        let (_, mut g_x) = Self::div_rem(self, f);

        while k > 0 {
            if k & 1 == 1 {
                let prod = &s * &g_x;
                let (_, rem) = Self::div_rem(&prod, f);
                s = rem;
            }

            let sq = &g_x * &g_x;
            let (_, rem) = Self::div_rem(&sq, f);
            g_x = rem;

            k >>= 1;
        }

        s
    }

    // Handbook of Applied Cryptography, Algorithm 4.69, Testing a polynomial for irreducibility (Ben-Or)
    // INPUT: a prime p and a monic polynomial f of degree m over Z_p[x]
    // OUTPUT: an answer to the question "is f irreducible over Z_p[x]?"
    pub fn is_irreducible(&self) -> bool {
        if self.is_zero() || self.deg() == 0 {
            return false;
        }

        let mut u = poly!(0, 1);
        let q = 1usize << (M as usize);

        for _ in 1..=(self.deg() / 2) {
            u = u.mod_pow(q, self);

            let mut u_minus_x = u.clone();
            if u_minus_x.coeffs.len() > 1 {
                u_minus_x.coeffs[1] = u_minus_x.coeffs[1] + GF::new(1);
            } else {
                u_minus_x.coeffs = vec![u.coeffs.get(0).copied().unwrap_or(GF::new(0)), GF::new(1)];
            }
            u_minus_x.clean();

            let d = Polynomial::gcd(&u_minus_x, self);

            if d.deg() > 0 {
                return false;
            }
        }

        true
    }
}

impl<const M: u8> Index<usize> for Polynomial<M> {
    type Output = GF<M>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coeffs[index]
    }
}

impl<const M: u8> Mul for &Polynomial<M> {
    type Output = Polynomial<M>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = vec![GF::new(0); self.coeffs.len() + rhs.coeffs.len() - 1];

        for i in 0..self.coeffs.len() {
            for j in 0..rhs.coeffs.len() {
                res[i + j] = res[i + j] + (self[i] * rhs[j]);
            }
        }

        let mut res_poly = Polynomial::new(res);
        res_poly.clean();

        res_poly
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf::GF;

    type TestGF = GF<13>;
    type TestPoly = Polynomial<13>;

    #[test]
    fn test_poly_macro_and_clean() {
        let mut p: TestPoly = poly![1, 2, 0, 0];
        p.clean();
        assert_eq!(p.deg(), 1);
        assert_eq!(p.coeffs.len(), 2);
    }

    #[test]
    fn test_poly_eval() {
        let p: TestPoly = poly![1, 1, 1]; // x^2 + x + 1
        let x_val = TestGF::new(2);
        let result = p.eval(x_val);
        assert_eq!(result, TestGF::new(7));
    }

    #[test]
    fn test_poly_div_rem() {
        let p: TestPoly = poly![1, 1, 1, 1]; // x^3 + x^2 + x + 1
        let q: TestPoly = poly![1, 1]; // x + 1
        let (div, rem) = TestPoly::div_rem(&p, &q);
        assert_eq!(div, poly![1, 0, 1]); // x^2 + 1
        assert!(rem.is_zero());
    }

    #[test]
    fn test_poly_gcd() {
        let p: TestPoly = poly![1, 0, 0, 1]; // x^3 + 1
        let q: TestPoly = poly![1, 0, 1]; // x^2 + 1
        let gcd = TestPoly::gcd(&p, &q);
        assert_eq!(gcd, poly![1, 1]); // x + 1
    }

    #[test]
    fn test_poly_mod_pow() {
        let p: TestPoly = poly![0, 1]; // x
        let q: TestPoly = poly![1, 1, 1]; // x^2 + x + 1
        let result = p.mod_pow(3, &q);
        assert_eq!(result, poly![1]); //  1
        assert_eq!(result.deg(), 0);
    }

    #[test]
    fn test_poly_is_irreducible() {
        let p: TestPoly = poly![1, 1, 0, 1]; // x^3 + x + 1
        assert!(p.is_irreducible());
        let q: TestPoly = poly![1, 0, 0, 1]; // x^3 + 1 = (x + 1)(x^2 + x + 1)
        assert!(!q.is_irreducible());
    }
}
