use aes::Aes256;
use aes::cipher::{BlockEncrypt, KeyInit, generic_array::GenericArray};
use rand::RngCore;

/// Byte-exact port of NIST PQC KAT reference `rng.c` (AES256_CTR_DRBG)
pub struct NistAesCtrDrbg {
    key: [u8; 32],
    v: [u8; 16],
}

impl NistAesCtrDrbg {
    /// entropy_input is the 48-byte `seed` field from the .rsp file.
    pub fn new(entropy_input: &[u8; 48]) -> Self {
        let mut drbg = NistAesCtrDrbg {
            key: [0u8; 32],
            v: [0u8; 16],
        };
        drbg.update(Some(entropy_input));
        drbg
    }

    fn aes256_ecb_encrypt_block(key: &[u8; 32], block_in: &[u8; 16]) -> [u8; 16] {
        let cipher = Aes256::new(GenericArray::from_slice(key));
        let mut block = *GenericArray::from_slice(block_in);
        cipher.encrypt_block(&mut block);
        block.into()
    }

    // V is treated as a 16-byte big-endian counter (increment from the last byte)
    fn increment_v(&mut self) {
        for j in (0..16).rev() {
            if self.v[j] == 0xff {
                self.v[j] = 0x00;
            } else {
                self.v[j] += 1;
                break;
            }
        }
    }

    fn update(&mut self, provided_data: Option<&[u8; 48]>) {
        let mut temp = [0u8; 48];
        for i in 0..3 {
            self.increment_v();
            let block = Self::aes256_ecb_encrypt_block(&self.key, &self.v);
            temp[16 * i..16 * i + 16].copy_from_slice(&block);
        }
        if let Some(pd) = provided_data {
            for i in 0..48 {
                temp[i] ^= pd[i];
            }
        }
        self.key.copy_from_slice(&temp[0..32]);
        self.v.copy_from_slice(&temp[32..48]);
    }

    pub fn randombytes(&mut self, out: &mut [u8]) {
        let xlen = out.len();
        let mut i = 0;
        while i < xlen {
            self.increment_v();
            let block = Self::aes256_ecb_encrypt_block(&self.key, &self.v);
            let remaining = xlen - i;
            if remaining > 16 {
                out[i..i + 16].copy_from_slice(&block);
                i += 16;
            } else {
                out[i..xlen].copy_from_slice(&block[..remaining]);
                i = xlen;
            }
        }
        self.update(None);
    }
}

impl RngCore for NistAesCtrDrbg {
    fn next_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill_bytes(&mut buf);
        u32::from_le_bytes(buf)
    }
    fn next_u64(&mut self) -> u64 {
        let mut buf = [0u8; 8];
        self.fill_bytes(&mut buf);
        u64::from_le_bytes(buf)
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.randombytes(dest);
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}
