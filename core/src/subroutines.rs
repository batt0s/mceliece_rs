use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::mceliece::{Ciphertext, PrivateKey, PublicKey, SessionKey, SysGF, SysPoly};
use crate::params::{K_U64, PARAMS, POLY_CAPACITY};

type MatGenRes = (Vec<u64>, Vec<SysGF>);

// McEliece Specification (Section 4.2) Matrix Generation for Goppa codes (Systematic form)
// INPUT: g - generator polynomial, alphas - roots of the generator polynomial
// OUTPUT: Some((T, )) if successful, None otherwise
// TODO: Semi-systematic matrix generation
pub fn matgen(g: &SysPoly, alphas: &[SysGF]) -> Option<MatGenRes> {
    let m = PARAMS.m as usize;
    let n = PARAMS.n;
    let t = PARAMS.t;
    let mt = m * t;

    let n_u64 = (n + 63) / 64;

    let mut h_hat = vec![0u64; mt * n_u64];

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
                let row = i * m + k;

                let block = j / 64;
                let bit_index = j % 64;

                h_hat[row * n_u64 + block] |= (bit as u64) << bit_index;
            }
            num = num * a_j;
        }
    }

    let t_matrix = reduce_to_systematic_form(&mut h_hat, mt, n)?;

    Some((t_matrix, alphas.to_vec()))
}

fn reduce_to_systematic_form(matrix: &mut Vec<u64>, rows: usize, cols: usize) -> Option<Vec<u64>> {
    let k_cols = cols - rows; // n - mt
    let n_u64 = (cols + 63) / 64;
    let k_u64 = (k_cols + 63) / 64;

    let mut is_invertible = Choice::from(1u8);

    // Make the matrix row-echelon and reduced form (Gauss-Jordan elimination)
    for i in 0..rows {
        let pivot_block = i / 64;
        let pivot_bit_idx = i % 64;

        let pivot_bit_i = (matrix[i * n_u64 + pivot_block] >> pivot_bit_idx) & 1u64;
        let mut need_swap = Choice::from((1u64 ^ pivot_bit_i) as u8);

        for j in (i + 1)..rows {
            let pivot_bit_j = (matrix[j * n_u64 + pivot_block] >> pivot_bit_idx) & 1u64;
            let bit_j_choice = Choice::from(pivot_bit_j as u8);

            let do_swap = need_swap & bit_j_choice;

            for c in 0..n_u64 {
                let mut val_i = matrix[i * n_u64 + c];
                let mut val_j = matrix[j * n_u64 + c];
                u64::conditional_swap(&mut val_i, &mut val_j, do_swap);

                matrix[i * n_u64 + c] = val_i;
                matrix[j * n_u64 + c] = val_j;
            }

            need_swap = need_swap & !do_swap;
        }

        let final_pivot = (matrix[i * n_u64 + pivot_block] >> pivot_bit_idx) & 1u64;
        is_invertible = is_invertible & Choice::from(final_pivot as u8);

        for j in 0..rows {
            if i == j {
                continue;
            }

            let bit_j = (matrix[j * n_u64 + pivot_block] >> pivot_bit_idx) & 1u64;
            let do_xor = Choice::from(bit_j as u8);

            for c in 0..n_u64 {
                let current = matrix[j * n_u64 + c];
                let xor_val = matrix[i * n_u64 + c];
                let new_val = current ^ xor_val;
                matrix[j * n_u64 + c] = u64::conditional_select(&current, &new_val, do_xor);
            }
        }
    }

    if is_invertible.unwrap_u8() == 0 {
        return None;
    }

    let mut t_matrix = vec![0u64; rows * k_u64];
    for i in 0..rows {
        for c in 0..k_cols {
            let source_col = rows + c;
            let bit = (matrix[i * n_u64 + (source_col / 64)] >> (source_col % 64)) & 1u64;
            t_matrix[i * k_u64 + (c / 64)] |= bit << (c % 64);
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
    let k_u64 = K_U64;

    let mut c_bits = vec![0u8; mt];

    // 1. C = e_1
    for i in 0..mt {
        c_bits[i] = e[i];
    }

    let mut e2_blocks = vec![0u64; k_u64];
    for j in 0..k {
        let bit = e[mt + j] as u64;
        e2_blocks[j / 64] |= bit << (j % 64);
    }

    // 2. C = C xor T_j for each j such that e_{mt+j} = 1 (C = C + T * e_2)
    for i in 0..mt {
        let mut dot_product = 0u64;
        for c in 0..k_u64 {
            dot_product ^= pk.T[i * k_u64 + c] & e2_blocks[c];
        }

        let parity = (dot_product.count_ones() % 2) as u8;
        c_bits[i] ^= parity;
    }

    c_bits
}

// McEliece Specification (Section 6.2) Representation of objects as byte strings
// Pack bits into bytes (Little-endian)
pub fn pack_bits(bits: &[u8]) -> Vec<u8> {
    let mut bytes = vec![0u8; (bits.len() + 7) / 8];
    for (i, &bit) in bits.iter().enumerate() {
        bytes[i / 8] |= bit << (i % 8);
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
pub fn decode(c: &Ciphertext, sk: &PrivateKey) -> (Vec<u8>, Choice) {
    let n = PARAMS.n;
    let t = PARAMS.t;
    let mt = (PARAMS.m as usize) * t;

    // Step 1: Extend C to v = (C, 0, 0, ..., 0) in F_2^n by appending k zeros
    let mut v = unpack_bits(&c, mt);
    v.resize(n, 0);

    // Step 2: Find unique c in F_2^n such that Hc = 0 and c has Hamming distance <= t from v.
    // Step 2.1: Compute Syndrome S(x) = sum_{j=0}^{n-1} v_j / (x - alpha_j) mod g(x)
    let (s_poly, s_valid) = compute_syndrome(&v, &sk);

    // Step 2.2: Find sigma(x) error locator polynomial using the Patterson algorithm
    let (sigma_poly, patterson_valid) = patterson_error_locator(&s_poly, &sk.g);

    // Step 2.3: Find error positions using the error locator polynomial sigma(x) (Chien Search)
    let (e, is_chien_valid) = chien_search(&sigma_poly, &sk.alphas, n);

    let mut weight = 0usize;
    for i in 0..n {
        weight += e[i] as usize;
    }

    // Step 4: If wt(e) = t and C = He, return e. Otherwise return None.
    // if weight != t {
    //     return None;
    // }
    // Make it constant-time
    let is_weight_correct: Choice = weight.ct_eq(&t);

    // Verify that C = He
    // if !verify_syndrome(&e, sk, &s_poly) {
    //     return None;
    // }
    // Make it constant-time
    let is_syndrome_correct: Choice = verify_syndrome(&e, sk, &s_poly);

    let is_valid: Choice =
        is_weight_correct & is_syndrome_correct & is_chien_valid & patterson_valid & s_valid;

    (e, is_valid)
}

fn compute_syndrome(v: &[u8], sk: &PrivateKey) -> (SysPoly, Choice) {
    let mut s_poly = SysPoly::zero();
    let mut all_valid = Choice::from(1u8);

    for j in 0..PARAMS.n {
        let mut denom = SysPoly::zero();
        denom.coeffs[0] = sk.alphas[j];
        denom.coeffs[1] = SysGF::new(1);

        let (inv, is_invertible) = denom.inv_mod(&sk.g, PARAMS.t);
        all_valid = all_valid & is_invertible;

        let mask = Choice::from(v[j]);

        for i in 0..POLY_CAPACITY {
            s_poly.coeffs[i] = SysGF::conditional_select(
                &s_poly.coeffs[i],
                &(s_poly.coeffs[i] + inv.coeffs[i]),
                mask,
            );
        }
    }
    (s_poly, all_valid)
}

fn patterson_error_locator(s_poly: &SysPoly, g: &SysPoly) -> (SysPoly, Choice) {
    let mut valid = Choice::from(1u8);

    // 1. T(x) = S(x)^-1 mod g(x)
    let (t_poly, t_valid) = s_poly.inv_mod(g, PARAMS.t);
    valid = valid & t_valid;

    // 2. R(x) = T(x) + x mod g(x)
    let mut r_poly = t_poly;
    r_poly.coeffs[1] = r_poly.coeffs[1] + SysGF::new(1);
    r_poly = r_poly.reduce(g, PARAMS.t);

    // 3. Q(x) = sqrt(R(x)) mod g(x)
    let mut g_even = SysPoly::zero();
    let mut g_odd = SysPoly::zero();
    for i in 0..=g.deg() {
        if i % 2 == 0 {
            g_even.coeffs[i / 2] = g.coeffs[i].sqrt();
        } else {
            g_odd.coeffs[i / 2] = g.coeffs[i].sqrt();
        }
    }

    let (g_odd_inv, g_odd_valid) = g_odd.inv_mod(g, PARAMS.t);
    valid = valid & g_odd_valid;
    let sqrt_x = (&g_even * &g_odd_inv).reduce(g, PARAMS.t);

    // Extract odd and even R(x)
    let mut r_even = SysPoly::zero();
    let mut r_odd = SysPoly::zero();
    for i in 0..=r_poly.deg() {
        if i % 2 == 0 {
            r_even.coeffs[i / 2] = r_poly.coeffs[i].sqrt();
        } else {
            r_odd.coeffs[i / 2] = r_poly.coeffs[i].sqrt();
        }
    }

    // Q(x) = r_even + sqrt_x * r_odd mod g(x)
    let prod = &sqrt_x * &r_odd;
    let mut q = SysPoly::zero();

    for i in 0..=prod.deg() {
        q.coeffs[i] = q.coeffs[i] + prod.coeffs[i];
    }
    for i in 0..=r_even.deg() {
        q.coeffs[i] = q.coeffs[i] + r_even.coeffs[i];
    }
    q = q.reduce(g, PARAMS.t);

    // Solve a(x) * Q(x) = b(x) mod g(x) with Constant Time Extended Euclidean Algorithm
    let (newr, newa) = ct_patterson_eea(g, &q);

    // 5. sigma(x) = a(x)^2 + x * b(x)^2 (a = newr, b = newa)
    let mut sigma = SysPoly::zero();
    let stop_deg = PARAMS.t / 2;

    for i in 0..=stop_deg {
        sigma.coeffs[2 * i] = newr.coeffs[i].sq();
        sigma.coeffs[2 * i + 1] = newa.coeffs[i].sq();
    }

    (sigma, valid)
}

// Returns (c_vec, is_valid) (for constant time concerns)
fn chien_search(sigma: &SysPoly, alphas: &[SysGF], n: usize) -> (Vec<u8>, Choice) {
    let mut c_vec = vec![0u8; n];
    let mut root_count: u16 = 0;
    for j in 0..n {
        let eval_res = sigma.eval(alphas[j]).0;
        let is_root: Choice = eval_res.ct_eq(&0);
        c_vec[j] = u8::conditional_select(&0, &1, is_root);
        root_count = u16::conditional_select(&root_count, &root_count.wrapping_add(1), is_root)
    }
    let is_valid = root_count.ct_eq(&(PARAMS.t as u16));
    (c_vec, is_valid)
}

fn verify_syndrome(e: &[u8], sk: &PrivateKey, s_poly: &SysPoly) -> Choice {
    let (e_syndrome, _e_valid) = compute_syndrome(e, sk);
    let mut is_equal = Choice::from(1u8);
    for i in 0..PARAMS.t {
        is_equal = is_equal & e_syndrome.coeffs[i].0.ct_eq(&s_poly.coeffs[i].0);
    }
    is_equal
}

// Constant time bitonic sort
pub fn ct_sort(arr: &mut [(u32, u32)]) {
    let n = arr.len();
    let mut k = 2;

    while k <= n {
        let mut j = k / 2;
        while j > 0 {
            for i in 0..n {
                let l = i ^ j;
                if l > i {
                    let dir = (i & k) == 0;

                    let a_val = arr[i].0;
                    let b_val = arr[l].0;

                    // is b_val greater than a_val?
                    let (_, borrow) = b_val.overflowing_sub(a_val);
                    let is_greater = Choice::from(borrow as u8);

                    let mut should_swap = is_greater;
                    if !dir {
                        should_swap = !should_swap;
                    }

                    let temp_val_i = u32::conditional_select(&arr[i].0, &arr[l].0, should_swap);
                    let temp_val_l = u32::conditional_select(&arr[l].0, &arr[i].0, should_swap);

                    let temp_idx_i = u32::conditional_select(&arr[i].1, &arr[l].1, should_swap);
                    let temp_idx_l = u32::conditional_select(&arr[l].1, &arr[i].1, should_swap);

                    arr[i] = (temp_val_i, temp_idx_i);
                    arr[l] = (temp_val_l, temp_idx_l);
                }
            }
            j /= 2;
        }
        k *= 2;
    }
}

// Constant Time Patterson EEA
pub fn ct_patterson_eea(g: &SysPoly, q: &SysPoly) -> (SysPoly, SysPoly) {
    let mut r0 = *g;
    let mut r1 = *q;

    let mut a0 = SysPoly::zero();
    let mut a1 = SysPoly::zero();
    a1.coeffs[0] = SysGF::new(1);

    // save snapshot value to find out when to stop
    let mut saved_a = SysPoly::zero();
    let mut saved_r = SysPoly::zero();
    let mut has_saved = Choice::from(0);

    let stop_deg = PARAMS.t / 2;

    // Maximum number of iterations is 2 * t
    for _ in 0..(2 * PARAMS.t) {
        let deg_r0 = r0.deg();
        let deg_r1 = r1.deg();

        // 1. Save snapshot: save for the first time if deg(r0) falls below stop_deg
        let is_r0_ready = Choice::from((deg_r0 <= stop_deg) as u8);
        let should_save = is_r0_ready & !has_saved;

        for i in 0..POLY_CAPACITY {
            saved_r.coeffs[i] =
                SysGF::conditional_select(&saved_r.coeffs[i], &r0.coeffs[i], should_save);
            saved_a.coeffs[i] =
                SysGF::conditional_select(&saved_a.coeffs[i], &a0.coeffs[i], should_save);
        }
        has_saved = has_saved | should_save;

        // 2. CT Swap: if deg(r0) < deg(r1), swap
        let is_r0_lesser = Choice::from((deg_r0 < deg_r1) as u8);
        let is_r1_zero = r1.is_zero();
        let do_swap = is_r0_lesser & !is_r1_zero;

        SysPoly::swap(&mut r0, &mut r1, do_swap);
        SysPoly::swap(&mut a0, &mut a1, do_swap);

        let deg_r0 = r0.deg();
        let deg_r1 = r1.deg();
        let diff = deg_r0 - deg_r1;

        // 3. lead_r0 * lead_r1^-1
        let lead_r0 = r0.lead();
        let lead_r1 = r1.lead();
        let safe_lead_r1 = SysGF::conditional_select(&lead_r1, &SysGF::new(1), is_r1_zero);
        let multiplier = lead_r0 * safe_lead_r1.inv();

        let mut shifted_r1 = r1;
        let mut shifted_a1 = a1;

        for i in 0..POLY_CAPACITY {
            shifted_r1.coeffs[i] = shifted_r1.coeffs[i] * multiplier;
            shifted_a1.coeffs[i] = shifted_a1.coeffs[i] * multiplier;
        }

        // 4. Shift the polynomial left
        for dst in 0..(PARAMS.t + 1) {
            let in_bounds = Choice::from((dst >= diff) as u8);
            let raw_src = dst.wrapping_sub(diff);
            let src = if raw_src < POLY_CAPACITY { raw_src } else { 0 };
            let r1_val = r1.coeffs[src] * multiplier;
            let a1_val = a1.coeffs[src] * multiplier;
            shifted_r1.coeffs[dst] = SysGF::conditional_select(&SysGF::new(0), &r1_val, in_bounds);
            shifted_a1.coeffs[dst] = SysGF::conditional_select(&SysGF::new(0), &a1_val, in_bounds);
        }

        // 5. r0 = r0 - shifted_r1
        let apply_sub = !is_r1_zero;
        for i in 0..POLY_CAPACITY {
            let new_r0_coeff = r0.coeffs[i] + shifted_r1.coeffs[i];
            let new_a0_coeff = a0.coeffs[i] + shifted_a1.coeffs[i];

            r0.coeffs[i] = SysGF::conditional_select(&r0.coeffs[i], &new_r0_coeff, apply_sub);
            a0.coeffs[i] = SysGF::conditional_select(&a0.coeffs[i], &new_a0_coeff, apply_sub);
        }
    }

    (saved_r, saved_a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf::GF;
    use crate::params::PK_SIZE;
    use crate::poly::Polynomial;

    #[test]
    fn test_reduce_to_systematic_form() {
        let raw_matrix = vec![
            vec![0, 1, 1, 1, 0, 0],
            vec![1, 0, 0, 0, 1, 0],
            vec![1, 1, 0, 0, 0, 1],
        ];

        let rows = 3;
        let cols = 6;

        let n_u64 = (cols + 63) / 64;

        let mut matrix = vec![0u64; rows * n_u64];
        for i in 0..rows {
            for j in 0..cols {
                if raw_matrix[i][j] == 1 {
                    matrix[i * n_u64 + (j / 64)] |= 1u64 << (j % 64);
                }
            }
        }

        let result = reduce_to_systematic_form(&mut matrix, rows, cols);

        assert!(result.is_some(), "Returned None");

        let t_matrix = result.unwrap();

        let k_cols = cols - rows;
        let k_u64 = (k_cols + 63) / 64;

        for i in 0..rows {
            for j in 0..rows {
                let bit = (matrix[i * n_u64 + (j / 64)] >> (j % 64)) & 1u64;
                if i == j {
                    assert_eq!(bit, 1, "Diagonal (pivot) {} {} must be 1", i, j);
                } else {
                    assert_eq!(bit, 0, "Non-diagonal {} {} must be 0", i, j);
                }
            }
        }

        assert_eq!(t_matrix.len(), rows * k_u64);
    }

    #[test]
    fn test_matgen_dimensions() {
        let t = PARAMS.t;
        let m = PARAMS.m as usize;
        let n = PARAMS.n;
        let pk_size = PK_SIZE;

        let g = Polynomial::from_slice(&[GF::new(1); PARAMS.t + 1]);

        let mut alphas = Vec::with_capacity(n);
        for i in 0..n {
            alphas.push(GF::new(i as u16));
        }

        let result = matgen(&g, &alphas);

        if let Some((t_matrix, out_alphas)) = result {
            let mt = m * t;
            let k = n - mt;

            assert_eq!(t_matrix.len(), mt);
            assert_eq!(t_matrix.len(), pk_size);
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
        let k_u64 = K_U64;

        let mut t_matrix = vec![0u64; mt * k_u64];
        for i in 0..mt {
            for j in 0..k {
                let bit = ((i + j) % 2) as u64;
                t_matrix[i * k_u64 + (j / 64)] = bit << (j % 64);
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
                c2[i],
                (pk.T[i * k_u64] & 1) as u8,
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
        let (decoded, is_valid) = decode(&ciphertext, &sk);
        assert_eq!(is_valid.unwrap_u8(), 1, "Decode should succeed");

        assert_eq!(decoded, e, "Decoded error vector must match original");
    }
}
