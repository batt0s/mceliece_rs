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

    // self mod f
    pub fn reduce(&self, f: &Self) -> Self {
        let mut rem = self.clone();
        rem.clean();
        let f_deg = f.deg();
        let f_lead_inv = f.coeffs.last().unwrap().inv();
        while !rem.is_zero() && rem.deg() >= f_deg {
            let deg_diff = rem.deg() - f_deg;
            let ratio = *rem.coeffs.last().unwrap() * f_lead_inv;
            for (i, &fc) in f.coeffs.iter().enumerate() {
                if fc.0 != 0 {
                    rem.coeffs[deg_diff + i] = rem.coeffs[deg_diff + i] + fc * ratio;
                }
            }
            rem.clean();
        }
        rem
    }

    // Divide-and-conquer product tree for polynomial multiplication, reduced mod f_y
    fn product_tree(factors: &[Self], f_y: &Self) -> Self {
        match factors.len() {
            0 => Self::new(vec![GF::new(1)]), // empty product = 1
            1 => factors[0].clone(),          // base case
            _ => {
                // split down the middle
                let mid = factors.len() / 2;
                let left = Self::product_tree(&factors[..mid], f_y);
                let right = Self::product_tree(&factors[mid..], f_y);
                // multiply and reduce mod f_y
                let prod = &left * &right;
                prod.reduce(f_y)
            }
        }
    }

    fn product_tree_ext(
        factors: &[Self], // the conjugates directly
        f_y: &Self,
    ) -> Vec<Self> {
        // outer poly coefficients
        match factors.len() {
            0 => vec![Self::new(vec![GF::new(1)])],
            1 => {
                // (X + conj) = [conj, 1]
                vec![factors[0].clone(), Self::new(vec![GF::new(1)])]
            }
            _ => {
                let mid = factors.len() / 2;
                let left = Self::product_tree_ext(&factors[..mid], f_y);
                let right = Self::product_tree_ext(&factors[mid..], f_y);

                // multiply left and right as outer polynomials
                // reuse your existing loop logic, just extracted here
                let zero = Self::new(vec![GF::new(0)]);
                let mut res = vec![zero; left.len() + right.len() - 1];
                for (i, ca) in left.iter().enumerate() {
                    for (j, cb) in right.iter().enumerate() {
                        let prod = ca * cb;
                        let rem = Self::reduce(&prod, f_y);
                        let r = &mut res[i + j];
                        let len = r.coeffs.len().max(rem.coeffs.len());
                        r.coeffs.resize(len, GF::new(0));
                        for (k, &c) in rem.coeffs.iter().enumerate() {
                            r.coeffs[k] = r.coeffs[k] + c;
                        }
                    }
                }
                res
            }
        }
    }

    // Square a polynomial in characteristic 2, reduced mod f.
    // (sum a_i * y^i)^2 = sum a_i^2 * y^(2i)  -- cross terms vanish
    fn frobenius_sq(p: &Self, f: &Self) -> Self {
        let deg = p.coeffs.len();
        let mut res = vec![GF::new(0); 2 * deg - 1];
        for (i, &c) in p.coeffs.iter().enumerate() {
            res[2 * i] = c.sq(); // squaring each GF element, coefficients go to even positions
        }
        let r = Polynomial::new(res);
        r.reduce(f)
    }

    // Apply Frobenius: p -> p^(2^M) mod f
    // Instead of mod_pow(2^M, f) which does 2^M iterations,
    // we do M squarings — from 8192 iterations down to 13.
    fn frobenius(&self, f: &Self) -> Self {
        let mut result = self.clone();
        for _ in 0..M {
            result = Polynomial::frobenius_sq(&result, f);
        }
        result
    }

    pub fn minpoly(&self, f_y: &Self) -> Self {
        // Collect conjugates via Frobenius: beta, beta^q, beta^(q^2), ...
        let mut conjugates: Vec<Self> = Vec::new();
        let mut current = self.clone();
        loop {
            if conjugates.iter().any(|c| c == &current) {
                break;
            }
            conjugates.push(current.clone());
            // current = current.mod_pow(q, f_y);
            current = current.frobenius(f_y);
        }

        // Multiply out (X - conj_0)(X - conj_1)...
        let acc = Polynomial::product_tree_ext(&conjugates, f_y);

        // At this point acc[i] should each be a degree-0 polynomial (a scalar in GF<M>)
        // because minpoly lands back in GF(2^M)[y] — extract those scalars
        let scalar_coeffs: Vec<GF<M>> = acc
            .iter()
            .map(|p| p.coeffs.get(0).copied().unwrap_or(GF::new(0)))
            .collect();

        let mut result = Polynomial::new(scalar_coeffs);
        result.clean();
        result
    }

    // Extended Euclidean algorithm for polynomial inversion modulo f
    pub fn inv_mod(&self, f: &Self) -> Option<Self> {
        let mut t = Self::new(vec![GF::new(0)]);
        let mut newt = Self::new(vec![GF::new(1)]);
        let mut r = f.clone();
        let mut newr = self.clone();
        while !newr.is_zero() {
            let (q, rem) = Self::div_rem(&r, &newr);
            r = newr;
            newr = rem;

            // next_t = t - q * newt
            let q_newt = &q * &newt;
            let len = std::cmp::max(t.coeffs.len(), q_newt.coeffs.len());
            let mut next_t_coeffs = vec![GF::new(0); len];
            for (i, c) in t.coeffs.iter().enumerate() {
                next_t_coeffs[i] = *c;
            }
            for (i, c) in q_newt.coeffs.iter().enumerate() {
                next_t_coeffs[i] = next_t_coeffs[i] + *c;
            }

            let mut next_t = Self::new(next_t_coeffs);
            next_t.clean();

            t = newt;
            newt = next_t;
        }

        if r.deg() > 0 {
            return None;
        }

        let scalar_inv = r.coeffs[0].inv();
        for c in &mut t.coeffs {
            *c = *c * scalar_inv;
        }
        Some(t)
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
    const T: usize = 96;

    fn f_y_460896() -> TestPoly {
        // F(y) = y^96 + y^10 + y^9 + y^6 + 1, z = GF<13>(2)
        let mut coeffs = vec![TestGF::new(0); T + 1];
        coeffs[0] = TestGF::new(1); // 1
        coeffs[6] = TestGF::new(1); // y^6
        coeffs[9] = TestGF::new(1); // y^9
        coeffs[10] = TestGF::new(1); // y^10
        coeffs[96] = TestGF::new(1); // y^96
        TestPoly::new(coeffs)
    }

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

    #[test]
    fn test_poly_minpoly() {
        let f_y = f_y_460896();

        // beta = 1 + y + y^2
        let mut beta_coeffs = vec![TestGF::new(0); 3];
        beta_coeffs[0] = TestGF::new(1);
        beta_coeffs[1] = TestGF::new(1);
        beta_coeffs[2] = TestGF::new(1);
        let beta = TestPoly::new(beta_coeffs);

        let g = beta.minpoly(&f_y);

        assert_eq!(g.coeffs.last().unwrap().0, 1, "minpoly must be monic");
        assert!(
            T % g.deg() == 0,
            "deg(g) = {} must divide T = {}",
            g.deg(),
            T
        );
        assert!(g.is_irreducible(), "minpoly must be irreducible");

        // g is USE 1: a real polynomial with GF13 scalar coefficients
        // beta is USE 2: an extension ring element represented as Polynomial<13>

        // g(beta) means: for each term g[i] * y^i in g,
        //   substitute beta for y -> g[i] * beta^i
        //   where g[i] is a GF13 scalar  (scales the extension element)
        //   and   beta^i is computed via repeated mul + div_rem mod f_y
        let mut result = TestPoly::new(vec![TestGF::new(0)]);
        let mut beta_pow = TestPoly::new(vec![TestGF::new(1)]); // beta^0 = 1

        for i in 0..=g.deg() {
            let scaled: Vec<TestGF> = beta_pow.coeffs.iter().map(|&c| c * g.coeffs[i]).collect();
            let scaled_poly = TestPoly::new(scaled);

            let len = result.coeffs.len().max(scaled_poly.coeffs.len());
            let mut res = vec![TestGF::new(0); len];
            for (j, c) in result.coeffs.iter().enumerate() {
                res[j] = res[j] + *c;
            }
            for (j, c) in scaled_poly.coeffs.iter().enumerate() {
                res[j] = res[j] + *c;
            }
            result = TestPoly::new(res);
            result.clean();

            let prod = &beta_pow * &beta;
            let (_, rem) = TestPoly::div_rem(&prod, &f_y);
            beta_pow = rem;
        }

        result.clean();
        assert!(
            result.is_zero(),
            "g(beta) must be 0 in GF(2^13), got {:?}",
            result
        );
    }

    #[test]
    fn test_poly_inv_mod() {
        let f = f_y_460896();

        // a(x) = x^2 + x + 1
        let a = poly![1, 1, 1];

        let inv_opt = a.inv_mod(&f);
        assert!(inv_opt.is_some(), "a(x) must have an inverse modulo f(x)");

        let inv = inv_opt.unwrap();

        // a(x) * a^-1(x)
        let prod = &a * &inv;

        let rem = TestPoly::reduce(&prod, &f);

        assert_eq!(rem.deg(), 0, "a(x) * a^-1(x) mod f(x) should have degree 0");
        assert_eq!(rem.coeffs[0].0, 1, "a(x) * a^-1(x) mod f(x) must be 1");
    }
}
