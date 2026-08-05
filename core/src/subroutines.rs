use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::mceliece::{Ciphertext, PrivateKey, PublicKey, SysGF, SysPoly};
use crate::params::{K_U64, PARAMS, POLY_CAPACITY};

type MatGenRes = (Vec<u64>, Vec<SysGF>, Vec<usize>);

/// McEliece Specification (Section 4.2) Matrix Generation for Goppa codes.
///
/// Supports systematic form (mu=nu=0) and semi-systematic form (mu,nu>0).
///
/// # Arguments
/// * `g` - Goppa polynomial
/// * `alphas` - Field elements α₀,…,α_{n-1}; permuted in-place for semi-systematic
///
/// # Returns
/// `(MatGenRes, Choice)` where:
/// * `T_matrix` - mt × k public-key matrix (undefined if `!is_valid`)
/// * `pivot_cols` - length mt; for i ∈ [mt-µ, mt) holds the original pivot column cᵢ
/// * `alphas` - Permuted field elements (for semi-systematic form)
/// * `is_valid` - `Choice(1)` iff reduction succeeded
///
/// # Constant-time
/// No. Contains an early return when `g_a_j == 0`, and the semi-systematic
/// alpha permutation branch is data-dependent.
pub fn matgen(g: &SysPoly, mut alphas: Vec<SysGF>) -> (MatGenRes, Choice) {
    let m = PARAMS.m as usize;
    let n = PARAMS.n;
    let t = PARAMS.t;
    let mt = m * t;
    let n_u64 = n.div_ceil(64);

    // Build h_hat matrix
    let mut h_hat = vec![0u64; mt * n_u64];
    for j in 0..n {
        let a_j = alphas[j];
        let g_a_j = g.eval(a_j);

        if g_a_j.0 == 0 {
            return (
                (vec![0u64; mt * K_U64], alphas, vec![0usize; mt]),
                Choice::from(0u8),
            );
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

    // Reduce to systematic form
    let (t_matrix, pivots, is_valid) =
        reduce_to_systematic_form(&mut h_hat, mt, n, PARAMS.mu, PARAMS.nu);

    // Permute alphas according to pivot_cols
    if PARAMS.semi_systematic {
        for (i, &pivot) in pivots
            .iter()
            .enumerate()
            .take(mt)
            .skip(mt.saturating_sub(PARAMS.mu))
        {
            alphas.swap(i, pivot);
        }
    }

    ((t_matrix, alphas, pivots), is_valid)
}

/// Pack the µ pivot-column offsets into the 8-byte `c` field of the private key.
///
/// Each original pivot column cᵢ (for i ∈ [mt-µ, mt)) is stored as the
/// *offset* `c_i - (mt - µ)`, which fits in ⌈log_2(ν)⌉ bits.
/// µ values are packed LSB-first into the 8-byte buffer.
///
/// # Constant-time
/// No. Contains an early return for `mu == 0` and uses data-dependent
/// iteration over the mu columns.
pub fn pack_col_perm(pivot_cols: &[usize]) -> [u8; 8] {
    let mu = PARAMS.mu;
    let mt = (PARAMS.m as usize) * PARAMS.t;

    if mu == 0 {
        return [255, 255, 255, 255, 0, 0, 0, 0];
    }

    let base = mt - mu;

    let mut pivots: u64 = 0;
    for i in 0..mu {
        let ci = pivot_cols[base + i];
        let offset = ci - base;
        pivots |= 1u64 << offset;
    }

    pivots.to_le_bytes()
}

fn reduce_to_systematic_form(
    matrix: &mut [u64],
    rows: usize,
    cols: usize,
    mu: usize,
    nu: usize,
) -> (Vec<u64>, Vec<usize>, Choice) {
    let k_cols = cols - rows; // n - mt
    let n_u64 = cols.div_ceil(64);
    let k_u64 = k_cols.div_ceil(64);

    let mut is_valid = Choice::from(1u8);
    let mut pivot_cols: Vec<usize> = (0..rows).collect();

    // Make the matrix row-echelon and reduced form (Gauss-Jordan elimination)
    for i in 0..rows {
        let search_start = i;
        let search_end = if i < rows.saturating_sub(mu) {
            i
        } else {
            let max_col = i.saturating_add(nu).saturating_sub(mu);
            max_col.min(cols.saturating_sub(1))
        };

        let mut found = Choice::from(0u8);
        let mut pivot_col = i;

        for col in search_start..=search_end {
            let c_block = col / 64;
            let c_bit = col % 64;

            let bit_i = (matrix[i * n_u64 + c_block] >> c_bit) & 1u64;
            let mut has_one = Choice::from(bit_i as u8);

            let do_scan = !found;
            for j in (i + 1)..rows {
                let bit_j = (matrix[j * n_u64 + c_block] >> c_bit) & 1u64;
                let bit_j_c = Choice::from(bit_j as u8);
                let do_swap = do_scan & bit_j_c & !has_one;

                for c in 0..n_u64 {
                    let mut vi = matrix[i * n_u64 + c];
                    let mut vj = matrix[j * n_u64 + c];
                    u64::conditional_swap(&mut vi, &mut vj, do_swap);
                    matrix[i * n_u64 + c] = vi;
                    matrix[j * n_u64 + c] = vj;
                }

                has_one |= do_scan & bit_j_c;
            }

            let this_is_pivot = has_one & !found;
            pivot_col = usize_cond_select(pivot_col, col, this_is_pivot);
            found |= this_is_pivot;
        }

        is_valid &= found;
        pivot_cols[i] = pivot_col;

        let eb = pivot_col / 64;
        let ebb = pivot_col % 64;
        for j in 0..rows {
            if i == j {
                continue;
            }
            let bit = (matrix[j * n_u64 + eb] >> ebb) & 1u64;
            let do_xor = Choice::from(bit as u8);
            for c in 0..n_u64 {
                let cur = matrix[j * n_u64 + c];
                let xor = matrix[i * n_u64 + c];
                matrix[j * n_u64 + c] = u64::conditional_select(&cur, &(cur ^ xor), do_xor);
            }
        }
    }

    for (i, &ci) in pivot_cols
        .iter()
        .enumerate()
        .take(rows)
        .skip(rows.saturating_sub(mu))
    {
        if ci == i {
            continue;
        }

        let i_block = i / 64;
        let c_block = ci / 64;
        let i_bit = i % 64;
        let c_bit = ci % 64;

        for row in 0..rows {
            let ptr = row * n_u64;
            let bit_i = (matrix[ptr + i_block] >> i_bit) & 1u64;
            let bit_c = (matrix[ptr + c_block] >> c_bit) & 1u64;

            matrix[ptr + i_block] &= !(1u64 << i_bit);
            matrix[ptr + c_block] &= !(1u64 << c_bit);
            matrix[ptr + i_block] |= bit_c << i_bit;
            matrix[ptr + c_block] |= bit_i << c_bit;
        }
    }

    let mut t_matrix = vec![0u64; rows * k_u64];
    for i in 0..rows {
        for c in 0..k_cols {
            let source_col = rows + c;
            let bit = (matrix[i * n_u64 + (source_col / 64)] >> (source_col % 64)) & 1u64;
            t_matrix[i * k_u64 + (c / 64)] |= bit << (c % 64);
        }
    }

    (t_matrix, pivot_cols, is_valid)
}

fn usize_cond_select(a: usize, b: usize, choice: Choice) -> usize {
    let mask = ((choice.unwrap_u8() as i8) as usize).wrapping_neg();
    a ^ ((a ^ b) & mask)
}

/// McEliece Specification (Section 4.3) Encoding Subroutine.
///
/// Encodes a weight-t column vector `e` in F_2^n using a public key `T`
/// into a ciphertext `C` in F_2^mt.
pub fn encode(e: &[u8], pk: &PublicKey) -> Vec<u8> {
    let mt = (PARAMS.m as usize) * PARAMS.t;
    let k = PARAMS.k;
    let k_u64 = K_U64;

    let mut c_bits = vec![0u8; mt];

    // 1. C = e_1
    c_bits[..mt].copy_from_slice(&e[..mt]);

    let mut e2_blocks = vec![0u64; k_u64];
    for j in 0..k {
        let bit = e[mt + j] as u64;
        e2_blocks[j / 64] |= bit << (j % 64);
    }

    // 2. C = C xor T_j for each j such that e_{mt+j} = 1 (C = C + T * e_2)
    for (i, c_bit) in c_bits.iter_mut().enumerate().take(mt) {
        let mut dot_product = 0u64;
        for (c, e2_block) in e2_blocks.iter().enumerate().take(k_u64) {
            dot_product ^= pk.T[i * k_u64 + c] & e2_block;
        }

        let parity = (dot_product.count_ones() % 2) as u8;
        *c_bit ^= parity;
    }

    c_bits
}

/// McEliece Specification (Section 6.2) Pack bits into bytes (Little-endian).
pub fn pack_bits(bits: &[u8]) -> Vec<u8> {
    let mut bytes = vec![0u8; bits.len().div_ceil(8)];
    for (i, &bit) in bits.iter().enumerate() {
        bytes[i / 8] |= bit << (i % 8);
    }
    bytes
}

/// Unpack bytes into bits (Little-endian).
pub fn unpack_bits(bytes: &[u8], num_bits: usize) -> Vec<u8> {
    let mut bits = Vec::with_capacity(num_bits);
    for i in 0..num_bits {
        bits.push((bytes[i / 8] >> (i % 8)) & 1);
    }
    bits
}

/// McEliece Specification (Section 4.4) Decoding subroutine.
///
/// Decodes `C` in F_2^mt to a word `e` of Hamming weight `wt(e) = t`
/// with `C = He` if such a word exists; otherwise returns failure.
///
/// # Constant-time
/// Yes. Returns `(Vec<u8>, Choice)` where the caller can use `Choice`
/// to mask the result. No early returns — all steps execute fully
/// regardless of intermediate validity.
pub fn decode(c: &Ciphertext, sk: &PrivateKey) -> (Vec<u8>, Choice) {
    let n = PARAMS.n;
    let t = PARAMS.t;
    let mt = (PARAMS.m as usize) * t;

    // Step 1: Extend C to v = (C, 0, 0, ..., 0) in F_2^n by appending k zeros
    let mut v = unpack_bits(c, mt);
    v.resize(n, 0);

    // Step 2: Find unique c in F_2^n such that Hc = 0 and c has Hamming distance <= t from v.
    // Step 2.1: Compute Syndrome S(x) = sum_{j=0}^{n-1} v_j / (x - alpha_j) mod g(x)
    let (s_poly, s_valid) = compute_syndrome(&v, sk);

    // Step 2.2: Find sigma(x) error locator polynomial using the Patterson algorithm
    let (sigma_poly, patterson_valid) = patterson_error_locator(&s_poly, &sk.g);

    // Step 2.3: Find error positions using the error locator polynomial sigma(x) (Chien Search)
    let (e, is_chien_valid) = chien_search(&sigma_poly, &sk.alphas, n);

    let mut weight = 0usize;
    for &bit in e.iter().take(n) {
        weight += bit as usize;
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

/// Computes the syndrome S(x) = Σ_j v_j · (x − α_j)⁻¹ mod g(x) of a word `v`.
///
/// # Fast inversion of linear factors
///
/// The textbook definition inverts (x − α_j) mod g(x) with a full polynomial
/// inversion for every position j, at Θ(t²) field multiplications each. For a
/// *linear* denominator there is a closed form. Let
///
///     h_j(x) = (g(x) − g(α_j)) / (x − α_j)
///
/// (well-defined because the numerator vanishes at x = α_j). Since
/// g(α_j) is a scalar, g(x) ≡ 0 (mod g(x)), and the field has characteristic
/// two (−1 = 1):
///
///     (x − α_j) · h_j(x) · g(α_j)⁻¹  =  g(x)·g(α_j)⁻¹ − 1  ≡  1  (mod g(x))
///
/// so (x − α_j)⁻¹ ≡ h_j(x) · g(α_j)⁻¹ (mod g(x)). Both g(α_j) and the
/// coefficients of h_j come out of one Horner pass (t multiplications), and
/// g(α_j)⁻¹ is a single constant-time scalar inversion (~2m multiplications).
///
/// # Constant-time
/// Yes. Every loop runs a fixed number of iterations (t per position) and
/// the g(α_j) = 0 case (non-invertible denominator) is handled with
/// conditional selection, mirroring `inv_mod`'s `is_invertible` flag.
fn compute_syndrome(v: &[u8], sk: &PrivateKey) -> (SysPoly, Choice) {
    let t = PARAMS.t;
    let g = &sk.g;
    let mut s_poly = SysPoly::zero();
    let mut all_valid = Choice::from(1u8);

    for (j, &v_j) in v.iter().enumerate().take(PARAMS.n) {
        let alpha = sk.alphas[j];

        // Single Horner pass computing h_j (degree t−1) and g(α_j):
        //   h_{t−1} = g_t,  h_{i−1} = g_i + α·h_i,  g(α) = g_0 + α·h_0
        // Requires g monic of degree t — guaranteed by keygen, which only
        // accepts polynomials with deg(g) == t from generate_irreducible.
        let mut h = SysPoly::zero();
        let mut acc = g.coeffs[t];
        h.coeffs[t - 1] = acc;
        for i in (1..t).rev() {
            acc = acc * alpha + g.coeffs[i];
            h.coeffs[i - 1] = acc;
        }
        let g_alpha = acc * alpha + g.coeffs[0];

        // If g(α_j) == 0 the inverse does not exist: mask the term to zero
        // and mark the syndrome invalid.
        let is_zero = g_alpha.ct_eq(&SysGF::new(0));
        let inv = SysGF::conditional_select(&g_alpha.inv(), &SysGF::new(0), is_zero);
        all_valid &= !is_zero;

        // S(x) += v_j · (x − α_j)⁻¹ mod g(x)
        let mask = Choice::from(v_j);
        for i in 0..t {
            s_poly.coeffs[i] = SysGF::conditional_select(
                &s_poly.coeffs[i],
                &(s_poly.coeffs[i] + h.coeffs[i] * inv),
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
    valid &= t_valid;

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
    valid &= g_odd_valid;
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
        is_equal &= e_syndrome.coeffs[i].0.ct_eq(&s_poly.coeffs[i].0);
    }
    is_equal
}

/// Constant-time bitonic sort for `(value, index)` pairs.
///
/// Sorts by value while preserving the associated index.
///
/// # Constant-time
/// Yes. The number of compare-exchange operations depends only on
/// `arr.len()`, not on the data values. All comparisons use
/// `overflowing_sub` and `conditional_select`.
pub fn ct_sort(arr: &mut [(u32, u32)]) {
    let n = arr.len();
    // Pad to the next power of two with sentinel values (u32::MAX, u32::MAX).
    // Bitonic sort requires the array length to be a power of two. For
    // arbitrary-length arrays (e.g., t=96 for mceliece460896), we pad with
    // sentinels that sort to the end, ensuring all real elements are correctly
    // compared and sorted.
    let next_pow2 = n.next_power_of_two();

    let mut padded: Vec<(u32, u32)> = Vec::with_capacity(next_pow2);
    padded.extend(arr.iter().copied());
    for _ in n..next_pow2 {
        padded.push((u32::MAX, u32::MAX));
    }

    let mut k = 2;
    while k <= next_pow2 {
        let mut j = k / 2;
        while j > 0 {
            for i in 0..next_pow2 {
                let l = i ^ j;
                if l > i {
                    let dir = (i & k) == 0;

                    let a_val = padded[i].0;
                    let b_val = padded[l].0;

                    // is b_val greater than a_val?
                    let (_, borrow) = b_val.overflowing_sub(a_val);
                    let is_greater = Choice::from(borrow as u8);

                    let mut should_swap = is_greater;
                    if !dir {
                        should_swap = !should_swap;
                    }

                    let temp_val_i =
                        u32::conditional_select(&padded[i].0, &padded[l].0, should_swap);
                    let temp_val_l =
                        u32::conditional_select(&padded[l].0, &padded[i].0, should_swap);

                    let temp_idx_i =
                        u32::conditional_select(&padded[i].1, &padded[l].1, should_swap);
                    let temp_idx_l =
                        u32::conditional_select(&padded[l].1, &padded[i].1, should_swap);

                    padded[i] = (temp_val_i, temp_idx_i);
                    padded[l] = (temp_val_l, temp_idx_l);
                }
            }
            j /= 2;
        }
        k *= 2;
    }

    // Copy the sorted real elements back
    arr[..n].copy_from_slice(&padded[..n]);
}

/// Constant Time Patterson Extended Euclidean Algorithm.
///
/// Solves `a(x) * Q(x) = b(x) (mod g)` by division steps.
/// Returns `(r, a)` where `r` and `a` are the saved polynomials
/// from the iteration where `deg(r0)` first drops below `t/2`.
///
/// # Constant-time
/// Yes. Runs for exactly `2 * t` iterations regardless of inputs.
/// All data-dependent decisions use conditional selection.
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
        has_saved |= should_save;

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

        let (result, _pivots, is_valid) = reduce_to_systematic_form(&mut matrix, rows, cols, 0, 0);
        assert!(is_valid.unwrap_u8() != 0, "Not valid");

        let t_matrix = result;

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

        let ((t_matrix, out_alphas, pivots), is_valid) = matgen(&g, alphas);
        if is_valid.unwrap_u8() != 0 {
            let mt = m * t;
            assert_eq!(t_matrix.len(), mt * K_U64);
            assert_eq!(t_matrix.len(), pk_size);
            assert_eq!(out_alphas.len(), n);
            assert_eq!(pivots.len(), mt)
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
