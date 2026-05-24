#[derive(Debug)]
pub struct McElieceParams {
    pub m: u8,
    pub n: usize,
    pub t: usize,
    pub q: usize,
    pub k: usize,
    pub f_z: u32,
    pub f_y_tpye: u8,
    pub l: usize,
    pub factors: &'static [usize],
}

#[cfg(feature = "mceliece348864")]
pub const PARAMS: McElieceParams = McElieceParams {
    m: 12,
    n: 3488,
    t: 64,
    q: 1usize << 12,
    k: 3488 - (12 * 64),
    f_z: 0x1009, // z^12 + z^3 + 1
    f_y_tpye: 0,
    l: 256,
    factors: &[2],
};
