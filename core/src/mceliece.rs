use crate::gf::GF;
use crate::params::{PARAMS, POLY_CAPACITY};
use crate::poly::Polynomial;
use crate::subroutines::{ct_sort, decode, encode, matgen, pack_bits, pack_col_perm, unpack_bits};
use rand::{Rng, RngCore};
use serde::{Deserialize, Serialize};
use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// The field GF(2^m) element type for the active parameter set.
pub type SysGF = GF<{ PARAMS.m }>;
/// The polynomial type over GF(2^m) for the active parameter set.
pub type SysPoly = Polynomial<{ PARAMS.m }, POLY_CAPACITY>;

/// Ciphertext produced by encapsulation.
pub type Ciphertext = Vec<u8>;
/// Session key produced by (de)capsulation — 256-bit (32-byte) output.
pub type SessionKey = [u8; 32];

/// Classic McEliece public key.
///
/// Contains the compressed public-key matrix T (mt × k bits, stored in
/// 64-bit words). Can be serialized via [`to_bytes`](PublicKey::to_bytes).
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
#[allow(non_snake_case)]
pub struct PublicKey {
    /// Public-key matrix T (mt rows of k bits each, packed into u64 words).
    pub T: Vec<u64>,
}

impl PublicKey {
    /// Classic McEliece spec (Section 6.2) canonical byte representation.
    /// Each row of T is packed to ceil(k/8) bytes, LSB-first.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mt = (PARAMS.m as usize) * PARAMS.t;
        let k = PARAMS.k;
        let k_u64 = crate::params::K_U64;
        let row_bytes = k.div_ceil(8);

        let mut out = vec![0u8; mt * row_bytes];
        for row in 0..mt {
            for bit_idx in 0..k {
                let word = self.T[row * k_u64 + bit_idx / 64];
                let bit = (word >> (bit_idx % 64)) & 1;
                out[row * row_bytes + bit_idx / 8] |= (bit as u8) << (bit_idx % 8);
            }
        }
        out
    }

    /// Deserialize a public key from its canonical byte representation.
    ///
    /// Parses the packed T matrix (mt rows of ceil(k/8) bytes each, LSB-first).
    /// Returns `None` if the input has the wrong length.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        let mt = (PARAMS.m as usize) * PARAMS.t;
        let k = PARAMS.k;
        let k_u64 = crate::params::K_U64;
        let row_bytes = k.div_ceil(8);

        let expected_len = mt * row_bytes;
        if bytes.len() != expected_len {
            return None;
        }

        let mut t_mat = vec![0u64; mt * k_u64];
        for row in 0..mt {
            for bit_idx in 0..k {
                let byte_val = bytes[row * row_bytes + bit_idx / 8];
                let bit = (byte_val >> (bit_idx % 8)) & 1;
                t_mat[row * k_u64 + bit_idx / 64] |= (bit as u64) << (bit_idx % 64);
            }
        }

        Some(PublicKey { T: t_mat })
    }
}

/// Classic McEliece private key.
///
/// Contains the seed `delta`, semi-systematic column pack `c`,
/// Goppa polynomial `g`, field ordering `alphas`, and random string `s`.
pub struct PrivateKey {
    /// Seed used to (re)generate this key.
    pub delta: [u8; 32],
    /// Semi-systematic column permutation info (8 bytes).
    pub c: [u8; 8],
    /// Goppa polynomial g(x) of degree t.
    pub g: SysPoly,
    /// Field ordering — the n field elements α₀,…,α_{n-1}.
    pub alphas: Vec<SysGF>,
    /// Random string s of length ceil(n/8) bytes.
    pub s: Vec<u8>,
}

impl PrivateKey {
    /// Partial Classic McEliece spec (Section 6.2) secret key encoding.
    ///
    /// Serializes the private key as: delta || c || g_coeffs || controlbits || s.
    pub fn to_bytes(&self) -> Vec<u8> {
        let t = PARAMS.t;
        let m = PARAMS.m;

        let mut out = Vec::new();
        out.extend_from_slice(&self.delta);
        out.extend_from_slice(&self.c);

        for i in 0..t {
            out.extend_from_slice(&self.g.coeffs[i].0.to_le_bytes());
        }

        // Field ordering, encoded per Benes-network control bits
        let pi: Vec<u32> = self.alphas.iter().map(|a| reverse_bits(a.0, m)).collect();
        let cb_bits = controlbits(&pi);
        out.extend_from_slice(&crate::subroutines::pack_bits(&cb_bits));

        out.extend_from_slice(&self.s);
        out
    }

    /// Deserialize a private key from its canonical byte representation.
    ///
    /// Parses the format: delta (32) || c (8) || g_coeffs (t×2) || controlbits || s.
    /// Returns `None` if the input has the wrong length.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        let m = PARAMS.m as usize;
        let t = PARAMS.t;
        let q = PARAMS.q;
        let n = PARAMS.n;

        let n_bytes = n.div_ceil(8);

        // Control bits size: C(q) = q/2 * (2*m - 1) raw bits for a Benes network
        // permuting q = 2^m elements.
        let cb_bits_count = (q / 2) * (2 * m - 1);
        let cb_bytes = cb_bits_count.div_ceil(8);

        // g coefficients: degree t, monic (leading coefficient = 1 is implicit)
        // The serialization stores coefficients 0..t-1 (t elements).
        let g_bytes = t * 2;

        let expected_len = 32 + 8 + g_bytes + cb_bytes + n_bytes;
        if bytes.len() != expected_len {
            return None;
        }

        let mut offset = 0;

        // delta: 32 bytes
        let mut delta = [0u8; 32];
        delta.copy_from_slice(&bytes[offset..offset + 32]);
        offset += 32;

        // c: 8 bytes (pivot-column bitmap, stored as-is)
        let mut c = [0u8; 8];
        c.copy_from_slice(&bytes[offset..offset + 8]);
        offset += 8;

        // g coefficients: t × 2 bytes each (u16 little-endian).
        // Only coefficients 0..t-1 are stored; the monic leading
        // coefficient (y^t) is implicitly 1.
        let mut g_coeffs = [SysGF::new(0); POLY_CAPACITY];
        for coeff in g_coeffs.iter_mut().take(t) {
            let val = u16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
            *coeff = SysGF::new(val);
            offset += 2;
        }
        g_coeffs[t] = SysGF::new(1); // monic leading coefficient
        let g = Polynomial::new(g_coeffs);

        // Control bits: unpack from bytes, decode to permutation, convert to alphas
        let cb_packed = &bytes[offset..offset + cb_bytes];
        let cb_bits = crate::subroutines::unpack_bits(cb_packed, cb_bits_count);
        let pi = controlbits_decode(&cb_bits, q);
        let mut alphas = Vec::with_capacity(q);
        for &p in &pi {
            let alpha_val = reverse_bits(p as u16, PARAMS.m) as u16;
            alphas.push(SysGF::new(alpha_val));
        }
        offset += cb_bytes;

        // s: n_bytes
        let s = bytes[offset..offset + n_bytes].to_vec();

        Some(PrivateKey {
            delta,
            c,
            g,
            alphas,
            s,
        })
    }
}

/// Reverses the order of the low `m` bits of `value`. Self-inverse.
fn reverse_bits(value: u16, m: u8) -> u32 {
    let mut result = 0u32;
    for j in 0..m {
        let bit = (value >> j) & 1;
        result |= (bit as u32) << (m - 1 - j);
    }
    result
}

/// composeinv(c, pi) = c . pi^-1, per the spec's Python `composeinv`.
fn composeinv(c: &[u32], pi: &[u32]) -> Vec<u32> {
    let mut pairs: Vec<(u32, u32)> = pi.iter().copied().zip(c.iter().copied()).collect();
    ct_sort(&mut pairs);
    pairs.into_iter().map(|(_, y)| y).collect()
}

/// Classic McEliece spec `controlbits` algorithm. (https://classic.mceliece.org/mceliece-sage-20221023/controlbits.py.html)
/// pi must be a permutation of {0, ..., n-1} with n = 2^m.
fn controlbits(pi: &[u32]) -> Vec<u8> {
    let n = pi.len();
    let mut m = 1u32;
    while (1usize << m) < n {
        m += 1;
    }
    assert_eq!(
        1usize << m,
        n,
        "controlbits input length must be a power of two"
    );

    if m == 1 {
        return vec![pi[0] as u8];
    }

    let idx: Vec<u32> = (0..n as u32).collect();

    let p0: Vec<u32> = idx.iter().map(|&x| pi[(x ^ 1) as usize]).collect();
    let q0: Vec<u32> = idx.iter().map(|&x| pi[x as usize] ^ 1).collect();

    let piinv = composeinv(&idx, pi);

    // p, q = composeinv(p,q), composeinv(q,p)   [python line 1]
    let p1 = composeinv(&p0, &q0);
    let q1 = composeinv(&q0, &p0);

    let mut c: Vec<u32> = idx.iter().map(|&x| x.min(p1[x as usize])).collect();

    // p, q = composeinv(p,q), composeinv(q,p)   [python line 2, reusing p1,q1]
    let mut p = composeinv(&p1, &q1);
    let mut q = composeinv(&q1, &p1);

    for _ in 1..(m - 1) {
        let cp = composeinv(&c, &q);
        let new_p = composeinv(&p, &q);
        let new_q = composeinv(&q, &p);
        p = new_p;
        q = new_q;
        c = idx
            .iter()
            .map(|&x| c[x as usize].min(cp[x as usize]))
            .collect();
    }

    let half = n / 2;
    let f: Vec<u8> = (0..half).map(|j| (c[2 * j] % 2) as u8).collect();
    let big_f: Vec<u32> = idx
        .iter()
        .map(|&x| x ^ (f[(x as usize) / 2] as u32))
        .collect();
    let f_pi = composeinv(&big_f, &piinv);

    let l: Vec<u8> = (0..half).map(|k| (f_pi[2 * k] % 2) as u8).collect();
    let big_l: Vec<u32> = idx
        .iter()
        .map(|&y| y ^ (l[(y as usize) / 2] as u32))
        .collect();

    let m_arr = composeinv(&f_pi, &big_l);

    let sub_m: [Vec<u32>; 2] = [
        (0..half).map(|j| m_arr[2 * j] / 2).collect(),
        (0..half).map(|j| m_arr[2 * j + 1] / 2).collect(),
    ];

    let subz0 = controlbits(&sub_m[0]);
    let subz1 = controlbits(&sub_m[1]);

    let mut z = Vec::with_capacity(subz0.len() + subz1.len());
    for i in 0..subz0.len() {
        z.push(subz0[i]);
        z.push(subz1[i]);
    }

    let mut result = f;
    result.extend(z);
    result.extend(l);
    result
}

/// Inverse of `controlbits`.
///
/// Given the raw control bits (each element is 0 or 1) and the permutation
/// size n (must be a power of two), reconstructs the permutation pi such
/// that `controlbits(pi) == bits`.
fn controlbits_decode(bits: &[u8], n: usize) -> Vec<u32> {
    let mut m = 1u32;
    while (1usize << m) < n {
        m += 1;
    }
    assert_eq!(
        1usize << m,
        n,
        "controlbits_decode input size must be a power of two"
    );

    if m == 1 {
        // n = 2: a single control bit.
        // bits[0] == 0 → pi = [0, 1]; bits[0] == 1 → pi = [1, 0]
        let b = bits[0] as u32;
        return vec![b, b ^ 1];
    }

    let half = n / 2;

    // Structure of the bit stream: f (n/2 bits) || z (recursive) || l (n/2 bits)
    let f = &bits[..half];
    let l = &bits[bits.len() - half..];
    let z = &bits[half..bits.len() - half];

    // De-interleave z: even indices → z0, odd indices → z1
    let sub_len = z.len() / 2;
    let mut z0 = Vec::with_capacity(sub_len);
    let mut z1 = Vec::with_capacity(sub_len);
    for i in 0..sub_len {
        z0.push(z[2 * i]);
        z1.push(z[2 * i + 1]);
    }

    // Recursively decode sub-permutations
    let sub_m0 = controlbits_decode(&z0, half);
    let sub_m1 = controlbits_decode(&z1, half);

    // Reconstruct m_arr from the decoded sub-permutations.
    // Encoding split: m_arr[2*j] / 2 → sub_m0[j], m_arr[2*j+1] / 2 → sub_m1[j]
    let mut m_arr = vec![0u32; n];
    for j in 0..half {
        m_arr[2 * j] = 2 * sub_m0[j];
        m_arr[2 * j + 1] = 2 * sub_m1[j] + 1;
    }

    // Reconstruct pi = big_F ∘ m_arr ∘ big_L.
    // Encoding: m_arr = big_F ∘ pi ∘ big_L   (composeinv(big_F, pi_inv) = big_F ∘ pi)
    // Both big_F and big_L are involutions (self-inverse), so:
    // pi = big_F ∘ m_arr ∘ big_L
    // where big_F[x] = x ^ f[x/2] and big_L[y] = y ^ l[y/2]
    let mut pi = vec![0u32; n];
    for i in 0..n {
        let after_l = (i as u32) ^ (l[i / 2] as u32);
        let after_m = m_arr[after_l as usize];
        let after_f = after_m ^ (f[(after_m / 2) as usize] as u32);
        pi[i] = after_f;
    }

    pi
}

/// Classic McEliece Specifications (Section 5.1) Irreducible-polynomial Generation.
///
/// Takes a string of `sigma_1 * t` input bits and outputs either `None` or
/// a monic irreducible degree-t polynomial g in F_q[x].
///
/// # Constant-time
/// No. Returns `None` (early return) when the bit length is wrong or when
/// `deg(g) != t`. This is a keygen function where the spec accepts
/// non-constant-time failure modes.
pub fn generate_irreducible(bits: &[u16]) -> Option<SysPoly> {
    let t = PARAMS.t;

    if bits.len() != PARAMS.t * (PARAMS.m as usize) {
        return None;
    }

    // Step 1 & 2: Build Beta_j scalars from bits, then build Beta polynomial
    let mut beta = SysPoly::zero();
    for j in 0..t {
        let mut beta_j = SysGF::new(0);
        for i in 0..(PARAMS.m as usize) {
            let bit_index = j * (PARAMS.m as usize) + i;
            if bits[bit_index] == 1 {
                beta_j.0 += SysGF::new(2).pow(i as u16).0;
            }
        }
        beta.coeffs[j] = beta_j;
    }

    // Step 3: Compute minimal polynomial of Beta in GF(2^M)
    let g = beta.minpoly(&PARAMS.f_y(), PARAMS.t);

    // Step 4: Return g if deg(g) == t
    if g.deg() == t { Some(g) } else { None }
}

/// Classic McEliece Specification (Section 5.2) Field-Ordering Generation.
///
/// Reads 32-bit values from `bytes`, sorts them using constant-time sort,
/// and generates the field ordering `alphas` via bit reversal.
///
/// # Constant-time
/// No. Returns `None` (early return) when the byte length is wrong or
/// when duplicate values are found. The `ct_sort` used internally is
/// constant-time, but the distinctness check uses an early-exit loop.
pub fn generate_field_ordering(bytes: &[u8]) -> Option<Vec<SysGF>> {
    let q = PARAMS.q;

    if bytes.len() != 4 * q {
        return None;
    }

    // Step 1: Read sigma_2 (32-bit) values
    let mut a: Vec<(u32, u32)> = Vec::with_capacity(q);
    for i in 0..q {
        let chunk = &bytes[4 * i..4 * (i + 1)];
        let a_i = u32::from_le_bytes(chunk.try_into().unwrap());
        a.push((a_i, i as u32));
    }

    // Step 2 & 3: Sort Lexicographically
    ct_sort(&mut a);

    // Check for distinction
    for i in 1..q {
        if a[i].0 == a[i - 1].0 {
            return None;
        }
    }

    // Step 4: Bit reversal from pi(i) indexes, generate alphas
    let mut alphas: Vec<SysGF> = Vec::with_capacity(q);
    for item in a.iter().take(q) {
        let pi_i = item.1;

        let mut alpha_val = 0u16;
        for j in 0..PARAMS.m {
            let bit = (pi_i >> j) & 1; // j. LSB
            alpha_val |= (bit as u16) << (PARAMS.m - 1 - j);
        }

        alphas.push(SysGF::new(alpha_val));
    }

    // Step 5: Return alphas
    Some(alphas)
}

/// Generates a key pair using system entropy.
///
/// Equivalent to `seeded_keygen(rand::thread_rng().fill(&mut seed))`.
///
/// # Constant-time
/// No. Key generation inherently uses rejection sampling loops and
/// may run for an unbounded number of iterations. Per the spec,
/// keygen timing leaks are an accepted trade-off.
pub fn keygen() -> (PublicKey, PrivateKey) {
    let mut seed = [0u8; 32];
    rand::thread_rng().fill(&mut seed);
    seeded_keygen(seed)
}

/// Deterministic key generation from a 32-byte seed.
///
/// Uses SHAKE-256 to derive randomness for the Goppa polynomial, field
/// ordering, and s-vector. Retries with a fresh seed if any step fails.
///
/// # Constant-time
/// No. Key generation uses rejection sampling loops and may run for
/// an unbounded number of iterations. Per the spec, keygen timing leaks
/// are an accepted trade-off.
pub fn seeded_keygen(mut seed: [u8; 32]) -> (PublicKey, PrivateKey) {
    let n = PARAMS.n;
    let t = PARAMS.t;
    let q = PARAMS.q;

    let n_bytes = n.div_ceil(8);
    let a_bytes = 4 * q;
    let g_bytes = 2 * t;

    loop {
        let mut hasher = Shake256::default();
        hasher.update(&[64u8]);
        hasher.update(&seed);
        let mut reader = hasher.finalize_xof();

        let mut s = vec![0u8; n_bytes];
        reader.read(&mut s);

        let mut a_buf = vec![0u8; a_bytes];
        reader.read(&mut a_buf);

        let mut g_buf = vec![0u8; g_bytes];
        reader.read(&mut g_buf);

        let mut next_seed = [0u8; 32];
        reader.read(&mut next_seed);

        let alphas = match generate_field_ordering(&a_buf) {
            Some(alphas) => alphas,
            None => {
                seed = next_seed;
                continue;
            }
        };

        let mut g_bits: Vec<u16> = Vec::with_capacity(t * (PARAMS.m as usize));
        for j in 0..t {
            // Convert 2 bytes to a 16-bit integer (little-endian)
            let block_16 = u16::from_le_bytes([g_buf[2 * j], g_buf[2 * j + 1]]);

            // Only take the first m bits and add to the bit vector
            for i in 0..(PARAMS.m as usize) {
                let bit = (block_16 >> i) & 1;
                g_bits.push(bit);
            }
        }

        let g = match generate_irreducible(&g_bits) {
            Some(g) => g,
            None => {
                seed = next_seed;
                continue;
            }
        };

        let ((t_matrix, alphas, pivots), is_valid) = matgen(&g, alphas);
        if is_valid.unwrap_u8() == 0 {
            seed = next_seed;
            continue;
        }

        let c = pack_col_perm(&pivots);

        let pk = PublicKey { T: t_matrix };
        let sk = PrivateKey {
            delta: seed,
            c,
            g,
            alphas,
            s,
        };
        return (pk, sk);
    }
}

/// Classic McEliece Specification (Section 5.4) Fixed-weight-vector generation.
///
/// Generates a random vector `e` in F_2^n with Hamming weight `t`.
/// Uses system entropy.
///
/// # Constant-time
/// No. Uses rejection sampling (potentially unbounded loop) to find
/// sufficient distinct indices less than n.
pub fn generate_fixed_weight() -> Vec<u8> {
    generate_fixed_weight_with_rng(&mut rand::thread_rng())
}

/// Fixed-weight-vector generation with a provided RNG.
///
/// Generates a vector `e` in F_2^n with Hamming weight `t`. The
/// inner loop uses constant-time distinctness checking via
/// `ct_sort`. However, the outer rejection loop is not constant-time.
///
/// # Constant-time
/// Partially. The distinctness check uses constant-time `ct_sort`
/// and `ct_eq`, but the outer rejection sampling loop may iterate
/// an unbounded number of times.
pub fn generate_fixed_weight_with_rng<R: RngCore>(rng: &mut R) -> Vec<u8> {
    let n = PARAMS.n;
    let t = PARAMS.t;
    let q = PARAMS.q;
    let m = PARAMS.m as usize;

    let mut tau = t;
    let mut bound = q;
    while bound > n {
        tau *= 2;
        bound /= 2;
    }
    let m_mask = (1u16 << m) - 1;

    loop {
        // Step 1: Generate sigma_1*tau uniform random bits b_0, b_1, ..., b_{sigma_1*tau-1}
        // Note: sigma_1 = 16 bit
        let mut buf = vec![0u8; 2 * tau];
        rng.fill_bytes(&mut buf);

        let mut a = Vec::with_capacity(t);

        // Step 2 & 3: Extract d_j (first m bits of every 16 bit chunk); if d_j < n, add to a. Repeat until a has t elements.
        for j in 0..tau {
            let d_j = u16::from_le_bytes([buf[2 * j], buf[2 * j + 1]]) & m_mask;
            if (d_j as usize) < n {
                a.push(d_j as usize);
                if a.len() == t {
                    break;
                }
            }
        }

        // If a does not have t elements, repeat from Step 1
        if a.len() < t {
            continue;
        }

        // Step 4: If not all distinct, repeat from Step 1
        let mut a_sorted: Vec<(u32, u32)> = a
            .iter()
            .enumerate()
            .map(|(i, &x)| (x as u32, i as u32))
            .collect();
        ct_sort(&mut a_sorted);
        let mut distinct = Choice::from(1u8);
        for i in 1..t {
            distinct &= !a_sorted[i].0.ct_eq(&a_sorted[i - 1].0);
        }
        if distinct.unwrap_u8() == 0 {
            continue;
        }

        // Step 5: Define e = (e_0, ..., e_{n-1}) in F_2^n as the weight-t vector such that e_a_i = 1 for each i.
        let mut e = vec![0u8; n];

        for secret_index in a {
            for (i, e_i) in e.iter_mut().enumerate() {
                let is_match = (i as u32).ct_eq(&(secret_index as u32));
                *e_i = u8::conditional_select(e_i, &1, is_match)
            }
        }

        // Step 6: Return e = (e_0, ..., e_{n-1})
        return e;
    }
}

/// McEliece Specification (Section 5.5) Encapsulation.
///
/// Takes a public key and outputs a ciphertext `C` and a 256-bit session key `K`.
/// Uses system entropy.
///
/// # Constant-time
/// Partially. Internally calls `generate_fixed_weight_with_rng` which uses
/// rejection sampling.
pub fn encapsulate(pk: &PublicKey) -> (Ciphertext, SessionKey) {
    encapsulate_with_rng(pk, &mut rand::thread_rng())
}

/// Encapsulation with a provided RNG.
///
/// # Arguments
/// * `pk` - The recipient's public key
/// * `rng` - A cryptographically secure RNG
///
/// # Returns
/// `(Ciphertext, SessionKey)` — the ciphertext and derived session key.
///
/// # Constant-time
/// Partially. Calls `generate_fixed_weight_with_rng` which uses rejection
/// sampling. The session key derivation via SHAKE-256 is constant-time.
pub fn encapsulate_with_rng<R: RngCore>(pk: &PublicKey, rng: &mut R) -> (Ciphertext, SessionKey) {
    // Step 1: Generate a random vector e with weight t.
    let e = generate_fixed_weight_with_rng(rng);

    // Step 2: Compute the ciphertext C = ENCODE(e, T)
    let c_bits = encode(&e, pk);

    // Step 3: Compute the session key K = H(1, e, C)
    let e_bytes = pack_bits(&e);
    let c_bytes = pack_bits(&c_bits);

    let mut hasher = Shake256::default();
    hasher.update(&[1u8]);
    hasher.update(&e_bytes);
    hasher.update(&c_bytes);

    let mut reader = hasher.finalize_xof();
    let mut k = [0u8; 32];
    reader.read(&mut k);

    // Step 4: Return (C, K)
    (c_bytes, k)
}

/// McEliece Specification (Section 5.6) Decapsulation.
///
/// Takes a ciphertext `C` and a private key, outputs a 256-bit session key.
///
/// If decoding fails, the session key is derived from the private key's `s`
/// string instead of the decoded error vector, preventing decryption failures
/// from leaking information.
///
/// # Constant-time
/// Yes. Regardless of whether decoding succeeds or fails, the function
/// produces a session key using the same code paths. The `is_valid`
/// `Choice` from `decode` is used to conditionally select between
/// `decoded_e` and `s_bits` using `u8::conditional_select`.
pub fn decapsulate(c: &Ciphertext, sk: &PrivateKey) -> SessionKey {
    // Step 1: Set b <- 1
    // Step 2: Extract s and Gamma' from private key
    let s = &sk.s;

    // Step 3: e <- DECODE(C, Gamma')
    let (decoded_e, is_valid) = decode(c, sk);

    let mut e = [0u8; PARAMS.n];

    // If decode fails (e_opt is None), e <- s and b <- 0
    let s_bits = unpack_bits(s, PARAMS.n);
    for i in 0..PARAMS.n {
        e[i] = u8::conditional_select(&s_bits[i], &decoded_e[i], is_valid);
    }
    let b = is_valid.unwrap_u8();

    // Step 4: Compute K = H(b, e, C)
    let e_bytes = pack_bits(&e);
    let mut hasher = Shake256::default();
    hasher.update(&[b]);
    hasher.update(&e_bytes);
    hasher.update(c);

    // Step 5: Return K
    let mut k = [0u8; 32];
    let mut reader = hasher.finalize_xof();
    reader.read(&mut k);

    k
}

#[cfg(test)]
mod tests {
    use crate::params::PK_SIZE;
    use proptest::prelude::*;
    use rand::SeedableRng;
    use std::sync::OnceLock;

    use super::*;

    /// Keypair shared by the roundtrip proptest. Keygen is expensive and
    /// seed-sensitive (seeded_keygen retries until the derived Goppa
    /// polynomial and field ordering are valid), so generate exactly once
    /// with a known-good seed: [7u8; 32] is also used by
    /// test_decapsulation_lifecycle.
    static KEYPAIR: OnceLock<(PublicKey, PrivateKey)> = OnceLock::new();

    #[test]
    fn test_generate_irreducible() {
        let bits: Vec<u16> = (0..PARAMS.t * 16).map(|i| (i % 2) as u16).collect();

        // May return None (wrong degree) — that's valid
        // If it returns Some, verify correctness
        if let Some(g) = generate_irreducible(&bits) {
            assert_eq!(g.deg(), PARAMS.t, "degree must be t");
            assert!(
                g.is_irreducible(g.deg()).unwrap_u8() == 1,
                "must be irreducible"
            );
        } else {
            println!("Returned None");
        }
    }

    #[test]
    fn test_generate_field_ordering() {
        let q = PARAMS.q;
        let mut bytes = Vec::with_capacity(q * 4);

        for i in 0..q {
            let a_i = (i as u32) * 17;
            bytes.extend_from_slice(&a_i.to_le_bytes());
        }

        let result = generate_field_ordering(&bytes);

        assert!(result.is_some(), "Returned None");
        let alphas = result.unwrap();
        assert_eq!(
            alphas.len(),
            q,
            "number of alphas must be equal to {}, got {:?}",
            q,
            alphas.len()
        );
    }

    #[test]
    fn test_keygen() {
        let (pk, sk) = keygen();

        let n = PARAMS.n;
        let t = PARAMS.t;
        let q = PARAMS.q;

        assert_eq!(
            pk.T.len(),
            PK_SIZE,
            "Size of T matrix should be {}",
            PK_SIZE
        );

        assert_eq!(sk.delta.len(), 32, "Delta (seed) should have 32 bytes");
        assert_eq!(
            sk.s.len(),
            (n + 7) / 8,
            "s vector should have ceil(n/8) bytes"
        );
        assert_eq!(sk.g.deg(), t, "g polynomial should have degree t");

        assert_eq!(sk.alphas.len(), q, "Alphas vector should have q elements");

        assert_eq!(sk.c.len(), 8, "c vector should have 8 bytes");
    }

    #[test]
    fn test_seeded_keygen_determinism() {
        let mut seed = [0u8; 32];
        for i in 0..32 {
            seed[i] = i as u8;
        }

        let (pk1, sk1) = seeded_keygen(seed);
        let (pk2, sk2) = seeded_keygen(seed);

        assert_eq!(pk1.T, pk2.T, "Public Key matrices should be the same");

        assert_eq!(sk1.delta, sk2.delta, "Delta (seed) should be the same");
        assert_eq!(sk1.c, sk2.c, "c vector should be constant");
        assert_eq!(sk1.s, sk2.s, "s vector should be the same");
        assert_eq!(sk1.alphas, sk2.alphas, "Alphas vector should be the same");

        assert_eq!(
            sk1.g.coeffs, sk2.g.coeffs,
            "Generated Irreducible polynomials should match"
        );
    }

    #[test]
    fn test_fixed_weight() {
        let n = PARAMS.n;
        let t = PARAMS.t;

        let e = generate_fixed_weight();

        assert_eq!(e.len(), n, "Vector length should be n");

        let mut weight = 0;
        for &bit in &e {
            assert!(bit == 0 || bit == 1, "Vector should only contain 0s and 1s");
            if bit == 1 {
                weight += 1;
            }
        }

        assert_eq!(weight, t, "Vector Hamming weight should be t");
    }

    #[test]
    fn test_encapsulate() {
        let (pk, _) = keygen();

        let (c_bytes, session_key) = encapsulate(&pk);

        let mt = (PARAMS.m as usize) * PARAMS.t;
        let expected_c_len = (mt + 7) / 8;

        assert_eq!(
            c_bytes.len(),
            expected_c_len,
            "Ciphertext byte length should be ceil(mt/8)"
        );
        assert_eq!(
            session_key.len(),
            32,
            "Session key should be 32 bytes (256 bits)"
        );

        assert_ne!(
            session_key, [0u8; 32],
            "Session key should not be all zeros"
        );
        assert_ne!(
            c_bytes,
            vec![0u8; expected_c_len],
            "Ciphertext should not be all zeros"
        );
    }

    #[test]
    fn test_decapsulation_lifecycle() {
        // Use seeded keygen for deterministic keys
        let seed = [7u8; 32];
        let (pk, sk) = seeded_keygen(seed);

        // Successful encapsulate/decapsulate roundtrip
        let (c_bytes, k_enc) = encapsulate(&pk);
        let k_dec = decapsulate(&c_bytes, &sk);
        assert_eq!(k_enc, k_dec, "Decapsulated key must match encapsulated key");

        // Tamper with ciphertext and ensure decapsulation produces a different key
        let mut c_tampered = c_bytes.clone();
        if !c_tampered.is_empty() {
            c_tampered[0] ^= 0xFF;
        }
        let k_dec_tampered = decapsulate(&c_tampered, &sk);
        assert_ne!(
            k_enc, k_dec_tampered,
            "Tampered ciphertext should not yield same session key"
        );
    }

    #[test]
    fn test_pk_roundtrip() {
        let (pk, _) = seeded_keygen([9u8; 32]);
        let bytes = pk.to_bytes();
        let pk2 = PublicKey::from_bytes(&bytes).expect("from_bytes should succeed");
        assert_eq!(
            pk.T, pk2.T,
            "Public key T matrices should match after round-trip"
        );
    }

    #[test]
    fn test_pk_from_bytes_wrong_len() {
        assert!(
            PublicKey::from_bytes(&[]).is_none(),
            "Empty bytes should return None"
        );
        assert!(
            PublicKey::from_bytes(&[0u8; 1]).is_none(),
            "Short bytes should return None"
        );
    }

    #[test]
    fn test_sk_roundtrip() {
        let (_, sk) = seeded_keygen([10u8; 32]);
        let bytes = sk.to_bytes();
        let sk2 = PrivateKey::from_bytes(&bytes).expect("from_bytes should succeed");

        assert_eq!(sk.delta, sk2.delta, "delta should match");
        assert_eq!(sk.c, sk2.c, "c should match");
        assert_eq!(sk.g.coeffs, sk2.g.coeffs, "g coefficients should match");
        assert_eq!(sk.alphas, sk2.alphas, "alphas should match");
        assert_eq!(sk.s, sk2.s, "s should match");
    }

    #[test]
    fn test_sk_from_bytes_wrong_len() {
        assert!(
            PrivateKey::from_bytes(&[]).is_none(),
            "Empty bytes should return None"
        );
    }

    #[test]
    fn test_full_serialization_roundtrip() {
        // Full life-cycle: keygen → serialize both keys → deserialize → encaps/decaps
        let (pk, sk) = seeded_keygen([11u8; 32]);

        let pk_bytes = pk.to_bytes();
        let sk_bytes = sk.to_bytes();

        let pk2 = PublicKey::from_bytes(&pk_bytes).expect("PK deserialize");
        let sk2 = PrivateKey::from_bytes(&sk_bytes).expect("SK deserialize");

        let (ct, k1) = encapsulate(&pk2);
        let k2 = decapsulate(&ct, &sk2);

        assert_eq!(
            k1, k2,
            "Session key should match after serialize/deserialize round-trip"
        );
    }

    // TODO: This takes around 5 minutes.
    proptest! {
        #![proptest_config(ProptestConfig::with_cases(5))]

        /// For the generated keypair, encapsulate(decapsulate(e)) must be the
        /// identity on the session key for every weight-t error vector.
        #[test]
        fn prop_encaps_decaps_roundtrip(seed in any::<[u8; 32]>()) {
            let (pk, sk) = seeded_keygen(seed);
            let (c, k_enc) = encapsulate(&pk);
            let k_dec = decapsulate(&c, &sk);
            prop_assert_eq!(
                k_enc, k_dec,
                "decapsulated key must match encapsulated key"
            );
        }
    }
}
