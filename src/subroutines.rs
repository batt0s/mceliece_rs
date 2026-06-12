use crate::mceliece::{Ciphertext, PrivateKey, PublicKey, SessionKey, SysGF, SysPoly};
use crate::params::PARAMS;

type MatGenRes = (Vec<Vec<u8>>, Vec<SysGF>);

// McEliece Specification (Section 4.2) Matrix Generation for Goppa codes (Systematic form)
// INPUT: g - generator polynomial, alphas - roots of the generator polynomial
// OUTPUT: Some((T, )) if successful, None otherwise
// TODO: Semi-systematic matrix generation
pub fn matgen(g: &SysPoly, alphas: &[SysGF]) -> Option<MatGenRes> {
    let m = PARAMS.m as usize;
    let n = PARAMS.n;
    let t = PARAMS.t;
    let mt = m * t;

    let mut h_hat = vec![vec![0u8; n]; mt];

    for j in 0..n {
        let a_j = alphas[j];
        let g_a_j = g.eval(a_j);

        if g_a_j.0 == 0 {
            return None;
        }

        let g_inv = g_a_j.inv();
        let mut num = SysGF::new(1);

        for i in 0..t {
            let h_ij = num * g_inv;
            for k in 0..m {
                let bit = ((h_ij.0 >> k) & 1) as u8;
                h_hat[i * m + k][j] = bit;
            }
            num = num * a_j;
        }
    }

    let t_matrix = reduce_to_systematic_form(&mut h_hat, mt, n)?;

    Some((t_matrix, alphas.to_vec()))
}

fn reduce_to_systematic_form(
    matrix: &mut Vec<Vec<u8>>,
    rows: usize,
    cols: usize,
) -> Option<Vec<Vec<u8>>> {
    let k_cols = cols - rows; // n - mt

    // Make the matrix row-echelon and reduced form (Gauss-Jordan elimination)
    for i in 0..rows {
        let mut pivot_row = i;
        let mut found = false;

        for j in i..rows {
            if matrix[j][i] == 1 {
                pivot_row = j;
                found = true;
                break;
            }
        }

        if !found {
            return None;
        }

        if pivot_row != i {
            matrix.swap(i, pivot_row);
        }

        let current_row = matrix[i].clone();

        for j in 0..rows {
            if j != i && matrix[j][i] == 1 {
                for c in i..cols {
                    matrix[j][c] ^= current_row[c];
                }
            }
        }
    }

    let mut t_matrix = vec![vec![0u8; k_cols]; rows];
    for i in 0..rows {
        for j in 0..k_cols {
            t_matrix[i][j] = matrix[i][rows + j];
        }
    }

    Some(t_matrix)
}

// McEliece Specification (Section 4.3) Encoding Subroutine
// INPUT: a weight-t column vector e in F_2^n and a public key T
// OUTPUT: a vector C in F_2^mt
pub fn encode(e: &[u8], pk: &PublicKey) -> Vec<u8> {
    let mt = (PARAMS.m as usize) * PARAMS.t;
    let k = PARAMS.k;

    let mut c_bits = vec![0u8; mt];

    // 1. C = e_1
    for i in 0..mt {
        c_bits[i] = e[i];
    }

    // 2. C = C xor T_j for each j such that e_{mt+j} = 1 (C = C + T * e_2)
    for j in 0..k {
        if e[mt + j] == 1 {
            for i in 0..mt {
                c_bits[i] ^= pk.T[i][j];
            }
        }
    }

    c_bits
}

// McEliece Specification (Section 6.2) Representation of objects as byte strings
// Pack bits into bytes (Little-endian)
pub fn pack_bits(bits: &[u8]) -> Vec<u8> {
    let mut bytes = vec![0u8; (bits.len() + 7) / 8];
    for (i, &bit) in bits.iter().enumerate() {
        if bit == 1 {
            bytes[i / 8] |= 1 << (i % 8);
        }
    }
    bytes
}

// Unpack bytes into bits (Little-endian)
pub fn unpack_bits(bytes: &[u8], num_bits: usize) -> Vec<u8> {
    let mut bits = Vec::with_capacity(num_bits);
    for i in 0..num_bits {
        bits.push((bytes[i / 8] >> (i % 8)) & 1);
    }
    bits
}

// McEliece Specification (Section 4.4) Decoding subroutine
// Decodes C in F_2^mt to a word e of Hamming weight wt(e) = t with C = He if such a word exists; otherwise it returns failure.
pub fn decode(c: &Ciphertext, sk: &PrivateKey) -> Option<Vec<u8>> {
    let n = PARAMS.n;
    let t = PARAMS.t;
    let mt = (PARAMS.m as usize) * t;

    // Step 1: Extend C to v = (C, 0, 0, ..., 0) in F_2^n by appending k zeros
    let mut v = unpack_bits(&c, mt);
    v.resize(n, 0);

    // Step 2: Find unique c in F_2^n such that Hc = 0 and c has Hamming distance <= t from v.
    // Step 2.1: Compute Syndrome S(x) = sum_{j=0}^{n-1} v_j / (x - alpha_j) mod g(x)
    let s_poly = compute_syndrome(&v, &sk);

    // Step 2.2: Find sigma(x) error locator polynomial using the Patterson algorithm
    let sigma_poly = patterson_error_locator(&s_poly, &sk.g)?;

    // Step 2.3: Find error positions using the error locator polynomial sigma(x) (Chien Search)
    let e_vec = chien_search(&sigma_poly, &sk.alphas, n);

    // If no error positions are found, return failure
    if e_vec.is_none() {
        return None;
    }
    let e = e_vec.unwrap();

    let mut weight = 0;
    for i in 0..n {
        if e[i] == 1 {
            weight += 1;
        }
    }

    // Step 4: If wt(e) = t and C = He, return e. Otherwise return None.
    if weight != t {
        return None;
    }

    // Verify that C = He
    if !verify_syndrome(&e, sk, &s_poly) {
        return None;
    }

    Some(e)
}

fn compute_syndrome(v: &[u8], sk: &PrivateKey) -> SysPoly {
    let mut s_poly = SysPoly::new(vec![SysGF::new(0)]);
    for j in 0..PARAMS.n {
        if v[j] == 1 {
            let denom = SysPoly::new(vec![sk.alphas[j], SysGF::new(1)]);
            if let Some(inv) = denom.inv_mod(&sk.g) {
                let len = std::cmp::max(s_poly.coeffs.len(), inv.coeffs.len());
                let mut next_s = vec![SysGF::new(0); len];
                for (i, c) in s_poly.coeffs.iter().enumerate() {
                    next_s[i] = next_s[i] + *c;
                }
                for (i, c) in inv.coeffs.iter().enumerate() {
                    next_s[i] = next_s[i] + *c;
                }
                s_poly = SysPoly::new(next_s);
                s_poly.clean();
            }
        }
    }
    s_poly
}

fn patterson_error_locator(s_poly: &SysPoly, g: &SysPoly) -> Option<SysPoly> {
    // 1. T(x) = S(x)^-1 mod g(x)
    let t_poly = s_poly.inv_mod(g)?;

    // 2. R(x) = T(x) + x mod g(x)
    let mut r_coeffs = t_poly.coeffs.clone();
    if r_coeffs.len() < 2 {
        r_coeffs.resize(2, SysGF::new(0));
    }
    r_coeffs[1] = r_coeffs[1] + SysGF::new(1);
    let mut r_poly = SysPoly::new(r_coeffs);
    r_poly.clean();
    r_poly = r_poly.reduce(g);

    // 3. Q(x) = sqrt(R(x)) mod g(x)
    let mut g_even = vec![SysGF::new(0); (g.coeffs.len() + 1) / 2];
    let mut g_odd = vec![SysGF::new(0); g.coeffs.len() / 2];
    for (i, &c) in g.coeffs.iter().enumerate() {
        if i % 2 == 0 {
            g_even[i / 2] = c.sqrt();
        } else {
            g_odd[i / 2] = c.sqrt();
        }
    }
    let g_odd_inv = SysPoly::new(g_odd).inv_mod(g)?;
    let sqrt_x = (&SysPoly::new(g_even) * &g_odd_inv).reduce(g);

    // Extract odd and even R(x)
    let mut r_even = vec![SysGF::new(0); (r_poly.coeffs.len() + 1) / 2];
    let mut r_odd = vec![SysGF::new(0); r_poly.coeffs.len() / 2];
    for (i, &c) in r_poly.coeffs.iter().enumerate() {
        if i % 2 == 0 {
            r_even[i / 2] = c.sqrt();
        } else {
            r_odd[i / 2] = c.sqrt();
        }
    }

    // Q(x) = r_even + sqrt_x * r_odd mod g(x)
    let prod = &sqrt_x * &SysPoly::new(r_odd);
    let len = std::cmp::max(prod.coeffs.len(), r_even.len());
    let mut q_coeffs = vec![SysGF::new(0); len];
    for (i, c) in prod.coeffs.iter().enumerate() {
        q_coeffs[i] = q_coeffs[i] + *c;
    }
    for (i, &c) in r_even.iter().enumerate() {
        q_coeffs[i] = q_coeffs[i] + c;
    }
    let mut q_poly: SysPoly = SysPoly::new(q_coeffs);
    q_poly.clean();
    q_poly = q_poly.reduce(g);

    // Solve a(x) * Q(x) = b(x) mod g(x) with Extended Euclidean Algorithm
    let mut a = SysPoly::new(vec![SysGF::new(0)]);
    let mut newa = SysPoly::new(vec![SysGF::new(1)]);
    let mut r = g.clone();
    let mut newr = q_poly;
    let stop_deg = PARAMS.t / 2;

    while newr.deg() > stop_deg {
        let (q_div, rem) = SysPoly::div_rem(&r, &newr);
        r = newr;
        newr = rem;

        let q_newa = &q_div * &newa;
        let len = std::cmp::max(q_newa.coeffs.len(), a.coeffs.len());
        let mut next_a = vec![SysGF::new(0); len];
        for (i, c) in a.coeffs.iter().enumerate() {
            next_a[i] = *c;
        }
        for (i, c) in q_newa.coeffs.iter().enumerate() {
            next_a[i] = next_a[i] + *c;
        }
        let mut next_a_poly = SysPoly::new(next_a);
        next_a_poly.clean();

        a = newa;
        newa = next_a_poly;
    }

    // 5. sigma(x) = a(x)^2 + x * b(x)^2 (a = newr, b = newa)
    let mut a_sq = vec![SysGF::new(0); 2 * newr.coeffs.len() - 1];
    for (i, &c) in newr.coeffs.iter().enumerate() {
        a_sq[2 * i] = c.sq();
    }

    let mut x_b_sq = vec![SysGF::new(0); 2 * newa.coeffs.len()];
    for (i, &c) in newa.coeffs.iter().enumerate() {
        x_b_sq[2 * i + 1] = c.sq();
    }

    let len = std::cmp::max(a_sq.len(), x_b_sq.len());
    let mut sigma_coeffs = vec![SysGF::new(0); len];
    for (i, &c) in a_sq.iter().enumerate() {
        sigma_coeffs[i] = sigma_coeffs[i] + c;
    }
    for (i, &c) in x_b_sq.iter().enumerate() {
        sigma_coeffs[i] = sigma_coeffs[i] + c;
    }

    let mut sigma = SysPoly::new(sigma_coeffs);
    sigma.clean();

    Some(sigma)
}

fn chien_search(sigma: &SysPoly, alphas: &[SysGF], n: usize) -> Option<Vec<u8>> {
    let mut c_vec = vec![0u8; n];
    let mut root_count = 0;
    for j in 0..n {
        if sigma.eval(alphas[j]).0 == 0 {
            c_vec[j] = 1;
            root_count += 1;
        }
    }
    if root_count != PARAMS.t {
        return None;
    }
    Some(c_vec)
}

fn verify_syndrome(e: &[u8], sk: &PrivateKey, s_poly: &SysPoly) -> bool {
    let e_syndrome = compute_syndrome(e, sk);
    e_syndrome.coeffs == s_poly.coeffs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf::GF;
    use crate::poly::Polynomial;

    #[test]
    fn test_reduce_to_systematic_form() {
        let mut matrix = vec![
            vec![0, 1, 1, 1, 0, 0],
            vec![1, 0, 0, 0, 1, 0],
            vec![1, 1, 0, 0, 0, 1],
        ];

        let rows = 3;
        let cols = 6;

        let result = reduce_to_systematic_form(&mut matrix, rows, cols);

        assert!(result.is_some(), "Returned None");

        let t_matrix = result.unwrap();

        for i in 0..rows {
            for j in 0..rows {
                if i == j {
                    assert_eq!(matrix[i][j], 1);
                } else {
                    assert_eq!(matrix[i][j], 0);
                }
            }
        }

        assert_eq!(t_matrix.len(), 3,);
        assert_eq!(t_matrix[0].len(), 3);
    }

    #[test]
    fn test_matgen_dimensions() {
        let t = PARAMS.t;
        let m = PARAMS.m as usize;
        let n = PARAMS.n;

        let g = Polynomial::new(vec![GF::new(1); t + 1]);

        let mut alphas = Vec::with_capacity(n);
        for i in 0..n {
            alphas.push(GF::new(i as u16));
        }

        let result = matgen(&g, &alphas);

        if let Some((t_matrix, out_alphas)) = result {
            let mt = m * t;
            let k = n - mt;

            assert_eq!(t_matrix.len(), mt);
            assert_eq!(t_matrix[0].len(), k);
            assert_eq!(out_alphas.len(), n);
        }
    }

    #[test]
    fn test_pack_bits() {
        let bits = vec![1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0];
        let packed = pack_bits(&bits);

        assert_eq!(packed.len(), 2, "11 bits should be packed into 2 bytes");
        assert_eq!(packed[0], 13, "0. byte has been corrupted");
        assert_eq!(packed[1], 3, "1. byte has been corrupted");

        let bits_8 = vec![1, 1, 1, 1, 1, 1, 1, 1];
        let packed_8 = pack_bits(&bits_8);
        assert_eq!(packed_8.len(), 1, "8 bits should be packed into 1 byte");
        assert_eq!(packed_8[0], 255, "packed byte should be 255 (all bits set)");
    }

    #[test]
    fn test_encode() {
        let mt = (PARAMS.m as usize) * PARAMS.t;
        let k = PARAMS.k;
        let n = PARAMS.n;

        let mut t_matrix = vec![vec![0u8; k]; mt];
        for i in 0..mt {
            for j in 0..k {
                t_matrix[i][j] = ((i + j) % 2) as u8;
            }
        }
        let pk = PublicKey {
            T: t_matrix.clone(),
        };

        let mut e1 = vec![0u8; n];
        e1[0] = 1;
        e1[3] = 1;
        let c1 = encode(&e1, &pk);
        assert_eq!(c1[0], 1, "First byte of C should be 1 (e_1 = [1, 0, 0, 1])");
        assert_eq!(c1[3], 1, "3. byte of C should be 1 (e_1 = [1, 0, 0, 1])");
        assert_eq!(
            c1.iter().sum::<u8>(),
            2,
            "Total number of set bits should be 2 (e_1 = [1, 0, 0, 1])"
        );

        let mut e2 = vec![0u8; n];
        e2[mt] = 1;
        let c2 = encode(&e2, &pk);
        for i in 0..mt {
            assert_eq!(
                c2[i], t_matrix[i][0],
                "C must match T matrix for e_2 = [0, 1, 0, 0]"
            );
        }
    }

    #[test]
    fn test_decode_roundtrip() {
        // Use seeded_keygen to obtain a matching keypair
        let seed = [42u8; 32];
        let (pk, sk) = crate::mceliece::seeded_keygen(seed);

        let n = PARAMS.n;
        let t = PARAMS.t;

        // Build an error vector e of weight t (first t positions set)
        let mut e = vec![0u8; n];
        for i in 0..t {
            e[i] = 1;
        }

        // Encode and pack into ciphertext
        let c_bits = encode(&e, &pk);
        let ciphertext = pack_bits(&c_bits);

        // Attempt to decode
        let decoded = decode(&ciphertext, &sk).expect("decode failed");

        assert_eq!(decoded, e, "Decoded error vector must match original");
    }
}
