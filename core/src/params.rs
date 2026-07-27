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
    pub semi_systematic: bool,
    pub mu: usize,
    pub nu: usize,
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
    semi_systematic: false,
    mu: 0,
    nu: 0,
};

#[cfg(feature = "mceliece348864f")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 12,
    n: 3488,
    t: 64,
    q: 1usize << 12,
    k: 3488 - (12 * 64),
    f_z: 0x1009, // z^12 + z^3 + 1
    f_y_coeffs: &F_Y_COEFFS_348864,
    l: 256,
    semi_systematic: true,
    mu: 32,
    nu: 64,
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
    semi_systematic: false,
    mu: 0,
    nu: 0,
};

#[cfg(feature = "mceliece460896f")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 4608,
    t: 96,
    q: 1usize << 13,
    k: 4608 - (13 * 96),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_460896,
    l: 256,
    semi_systematic: true,
    mu: 32,
    nu: 64,
};

#[cfg(feature = "mceliece6688128")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 6688,
    t: 128,
    q: 1usize << 13,
    k: 6688 - (13 * 128),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_6688128,
    l: 256,
    semi_systematic: false,
    mu: 0,
    nu: 0,
};

#[cfg(feature = "mceliece6688128f")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 6688,
    t: 128,
    q: 1usize << 13,
    k: 6688 - (13 * 128),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_6688128,
    l: 256,
    semi_systematic: true,
    mu: 32,
    nu: 64,
};

#[cfg(feature = "mceliece6960119")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 6960,
    t: 119,
    q: 1usize << 13,
    k: 6960 - (13 * 119),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_6960119,
    l: 256,
    semi_systematic: false,
    mu: 0,
    nu: 0,
};

#[cfg(feature = "mceliece6960119f")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 6960,
    t: 119,
    q: 1usize << 13,
    k: 6960 - (13 * 119),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_6960119,
    l: 256,
    semi_systematic: true,
    mu: 32,
    nu: 64,
};

#[cfg(feature = "mceliece8192128")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 8192,
    t: 128,
    q: 1usize << 13,
    k: 8192 - (13 * 128),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_8192128,
    l: 256,
    semi_systematic: false,
    mu: 0,
    nu: 0,
};

#[cfg(feature = "mceliece8192128f")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 13,
    n: 8192,
    t: 128,
    q: 1usize << 13,
    k: 8192 - (13 * 128),
    f_z: 0x201B, // z^13 + z^4 + z^3 + z + 1
    f_y_coeffs: &F_Y_COEFFS_8192128,
    l: 256,
    semi_systematic: true,
    mu: 32,
    nu: 64,
};

#[cfg(any(feature = "mceliece348864", feature = "mceliece348864f"))]
static F_Y_COEFFS_348864: [u16; 65] = {
    let mut arr = [0u16; 65];
    arr[0] = 2; // z (generator of GF(2^12))
    arr[1] = 1; // y
    arr[3] = 1; // y^3
    arr[64] = 1; // y^64
    arr
};

#[cfg(any(feature = "mceliece460896", feature = "mceliece460896f"))]
static F_Y_COEFFS_460896: [u16; 97] = {
    let mut arr = [0u16; 97];
    arr[0] = 1;
    arr[6] = 1;
    arr[9] = 1;
    arr[10] = 1;
    arr[96] = 1;
    arr
};

#[cfg(any(feature = "mceliece6688128", feature = "mceliece6688128f"))]
static F_Y_COEFFS_6688128: [u16; 129] = {
    let mut arr = [0u16; 129];
    arr[0] = 1;
    arr[1] = 1;
    arr[2] = 1;
    arr[7] = 1;
    arr[128] = 1;
    arr
};

#[cfg(any(feature = "mceliece6960119", feature = "mceliece6960119f"))]
static F_Y_COEFFS_6960119: [u16; 120] = {
    let mut arr = [0u16; 120];
    arr[0] = 1;
    arr[8] = 1;
    arr[119] = 1;
    arr
};

#[cfg(any(feature = "mceliece8192128", feature = "mceliece8192128f"))]
static F_Y_COEFFS_8192128: [u16; 129] = {
    let mut arr = [0u16; 129];
    arr[0] = 1;
    arr[1] = 1;
    arr[2] = 1;
    arr[7] = 1;
    arr[128] = 1;
    arr
};

#[cfg(any(feature = "mceliece348864", feature = "mceliece348864f"))]
pub const POLY_CAPACITY: usize = 128;

#[cfg(any(feature = "mceliece460896", feature = "mceliece460896f"))]
pub const POLY_CAPACITY: usize = 256;

#[cfg(any(feature = "mceliece6688128", feature = "mceliece6688128f"))]
pub const POLY_CAPACITY: usize = 256;

#[cfg(any(feature = "mceliece6960119", feature = "mceliece6960119f"))]
pub const POLY_CAPACITY: usize = 256;

#[cfg(any(feature = "mceliece8192128", feature = "mceliece8192128f"))]
pub const POLY_CAPACITY: usize = 256;

pub const MT: usize = (PARAMS.m as usize) * PARAMS.t;
pub const K_U64: usize = PARAMS.k.div_ceil(64);
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
