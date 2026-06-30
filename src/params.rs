use crate::gf::GF;
use crate::poly::Polynomial;

#[derive(Debug)]
pub struct McElieceParams {
    pub m: u8,
    pub n: usize,
    pub t: usize,
    pub q: usize,
    pub k: usize,
    pub f_z: u32,
    pub f_y_coeffs: &'static [u16],
    pub l: usize,
}

#[cfg(feature = "mceliece348864")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 12,
    n: 3488,
    t: 64,
    q: 1usize << 12,
    k: 3488 - (12 * 64),
    f_z: 0x1009, // z^12 + z^3 + 1
    f_y_coeffs: &F_Y_COEFFS_348864,
    l: 256,
};

static F_Y_COEFFS_348864: [u16; 65] = {
    let mut arr = [0u16; 65];
    arr[0] = 2; // z (generator of GF(2^12))
    arr[1] = 1; // y
    arr[3] = 1; // y^3
    arr[64] = 1; // y^64
    arr
};

#[cfg(feature = "mceliece460896")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 4608,
    t: 96,
    q: 1usize << 13,
    k: 4608 - (13 * 96),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_460896,
    l: 256,
};

static F_Y_COEFFS_460896: [u16; 97] = {
    let mut arr = [0u16; 97];
    arr[0] = 1;
    arr[6] = 1;
    arr[9] = 1;
    arr[10] = 1;
    arr[96] = 1;
    arr
};

pub const POLY_CAPACITY: usize = 256;
pub const MT: usize = (PARAMS.m as usize) * PARAMS.t;
pub const K_U64: usize = (PARAMS.k + 63) / 64;
pub const PK_SIZE: usize = MT * K_U64;

impl McElieceParams {
    pub fn f_y<const M: u8>(&self) -> Polynomial<M, POLY_CAPACITY> {
        let mut poly = Polynomial::zero();
        let len = self.f_y_coeffs.len().min(POLY_CAPACITY);
        for i in 0..len {
            poly.coeffs[i] = GF::<M>::new(self.f_y_coeffs[i]);
        }
        poly
    }
}
