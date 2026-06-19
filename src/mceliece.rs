use crate::gf::GF;
use crate::params::{PARAMS, POLY_CAPACITY};
use crate::poly::Polynomial;
use crate::subroutines::{decode, encode, matgen, pack_bits, unpack_bits};
use rand::Rng;
use serde::{Deserialize, Serialize};
use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};

pub type SysGF = GF<{ PARAMS.m }>;
pub type SysPoly = Polynomial<{ PARAMS.m }, POLY_CAPACITY>;

pub type Ciphertext = Vec<u8>;
pub type SessionKey = [u8; 32];

#[derive(Serialize, Deserialize)]
pub struct PublicKey {
    pub T: Vec<Vec<u8>>,
}

pub struct PrivateKey {
    pub delta: [u8; 32],
    pub c: [u8; 8],
    pub g: SysPoly,
    pub alphas: Vec<SysGF>,
    pub s: Vec<u8>,
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
                beta_j.0 = beta_j.0 + SysGF::new(2).pow(i as u16).0;
            }
        }
        beta.coeffs[j] = beta_j;
    }

    // Step 3: Compute minimal polynomial of Beta in GF(2^M)
    let g = beta.minpoly(&PARAMS.f_y());

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
    let mut a: Vec<(u32, usize)> = Vec::with_capacity(q);
    for i in 0..q {
        let chunk = &bytes[4 * i..4 * (i + 1)];
        let a_i = u32::from_le_bytes(chunk.try_into().unwrap());
        a.push((a_i, i));
    }

    // Step 2 & 3: Sort Lexicographically
    // TODO: Constant time sorting
    a.sort_unstable_by_key(|&(val, _)| val);

    // Check for distinction
    for i in 1..q {
        if a[i].0 == a[i - 1].0 {
            return None;
        }
    }

    // Step 4: Bit reversal from pi(i) indexes, generate alphas
    let mut alphas: Vec<SysGF> = Vec::with_capacity(q);
    for i in 0..q {
        let pi_i = a[i].1 as u32;

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

    let n_bytes = (n + 7) / 8;
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

        let (t_matrix, _) = match matgen(&g, &alphas) {
            Some(res) => res,
            None => {
                seed = next_seed;
                continue;
            }
        };

        let c: [u8; 8] = [255, 255, 255, 255, 0, 0, 0, 0];

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
    let n = PARAMS.n;
    let t = PARAMS.t;
    let q = PARAMS.q;
    let m = PARAMS.m as usize;

    // Determine tau based on n and t
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
        rand::thread_rng().fill(&mut buf[..]);

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
        let mut a_sorted = a.clone();
        a_sorted.sort_unstable();
        let mut distinct = true;
        for i in 1..t {
            if a_sorted[i] == a_sorted[i - 1] {
                distinct = false;
                break;
            }
        }
        if !distinct {
            continue;
        }

        // Step 5: Define e = (e_0, ..., e_{n-1}) in F_2^n as the weight-t vector such that e_a_i = 1 for each i.
        let mut e = vec![0u8; n];
        for i in a {
            e[i] = 1;
        }

        // Step 6: Return e = (e_0, ..., e_{n-1})
        return e;
    }
}

// McEliece Specification (Section 5.5) Encapsulation
// Takes a public key T. It outputs a ciphertext C and a session key K.
pub fn encapsulate(pk: &PublicKey) -> (Ciphertext, SessionKey) {
    // Step 1: Generate a random vector e with weight t.
    let e = generate_fixed_weight();

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
    let mut b = 1u8;

    // Step 2: Extract s and Gamma' from private key
    let s = &sk.s;

    // Step 3: e <- DECODE(C, Gamma')
    let e_opt = decode(c, &sk);

    // If decode fails (e_opt is None), e <- s and b <- 0
    let e = match e_opt {
        Some(e) => e,
        None => {
            b = 0;
            unpack_bits(s, PARAMS.n)
        }
    };

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
    use super::*;

    #[test]
    fn test_generate_irreducible() {
        let bits: Vec<u16> = (0..PARAMS.t * 16).map(|i| (i % 2) as u16).collect();

        // May return None (wrong degree) — that's valid
        // If it returns Some, verify correctness
        if let Some(g) = generate_irreducible(&bits) {
            assert_eq!(g.deg(), PARAMS.t, "degree must be t");
            assert!(g.is_irreducible(), "must be irreducible");
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
        let k = PARAMS.k; // n - mt
        let mt = m * t;

        assert_eq!(pk.T.len(), mt, "T matrix should have {} rows", mt);
        assert_eq!(pk.T[0].len(), k, "T matrix should have {} columns", k);

        assert_eq!(sk.delta.len(), 32, "Delta (seed) should have 32 bytes");
        assert_eq!(
            sk.s.len(),
            (n + 7) / 8,
            "s vector should have ceil(n/8) bytes"
        );
        assert_eq!(sk.g.deg(), t, "g polynomial should have degree t");

        assert_eq!(sk.alphas.len(), q, "Alphas vector should have q elements");

        assert_eq!(
            sk.c,
            [255, 255, 255, 255, 0, 0, 0, 0],
            "c vector should be constant"
        );
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
