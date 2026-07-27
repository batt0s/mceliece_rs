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

pub type SysGF = GF<{ PARAMS.m }>;
pub type SysPoly = Polynomial<{ PARAMS.m }, POLY_CAPACITY>;

pub type Ciphertext = Vec<u8>;
pub type SessionKey = [u8; 32];

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
#[allow(non_snake_case)]
pub struct PublicKey {
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
}

pub struct PrivateKey {
    pub delta: [u8; 32],
    pub c: [u8; 8],
    pub g: SysPoly,
    pub alphas: Vec<SysGF>,
    pub s: Vec<u8>,
}

impl PrivateKey {
    /// Partial Classic McEliece spec (Section 6.2) secret key encoding.
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

// Classic McEliece Specifications (Section 5.1) Irreducible-polynomial Generation
// The following algorithm Irreducible takes a string of sigma_1*t input bits d_0 , d_1 , . . . , d_{sigma_1 t−1} . It
// outputs either ⊥ or a monic irreducible degree-t polynomial g in F_q[x].
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

// Classic McEliece Specification (Section 5.2) Field-Ordering Generation
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

pub fn keygen() -> (PublicKey, PrivateKey) {
    let mut seed = [0u8; 32];
    rand::thread_rng().fill(&mut seed);
    seeded_keygen(seed)
}

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

// Classic McEliece Specification (Section 5.4) Fixed-weight-vector generation
// The following randomized algorithm takes no input.
// It outputs a vector e in F_2^n such that the Hamming weight of e is t.
// The algorithm uses a precomputed integer tau >= t. The integer tau is defined as
// if n = q; as 2t if q/2 <= n < q; as 4t if q/4 <= n < q/2; etc.
// All of the selected parameter sets have q/2 <= n < q, so tau in {t, 2t}.
pub fn generate_fixed_weight() -> Vec<u8> {
    generate_fixed_weight_with_rng(&mut rand::thread_rng())
}

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

// McEliece Specification (Section 5.5) Encapsulation
// Takes a public key T. It outputs a ciphertext C and a session key K.
pub fn encapsulate(pk: &PublicKey) -> (Ciphertext, SessionKey) {
    encapsulate_with_rng(pk, &mut rand::thread_rng())
}

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

// McEliece Specification (Section 5.6) Decapsulation
// Takes as input a ciphertext C and a private key, outputs a session key.
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

    use super::*;

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

        let m = PARAMS.m as usize;
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
}
