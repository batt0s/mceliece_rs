use mceliece::gf::GF;
use mceliece::params::PARAMS;
use mceliece::poly::Polynomial;

type SysGF = GF<{ PARAMS.m }>;
type SysPoly = Polynomial<{ PARAMS.m }>;

fn main() {
    println!("Selected parameters: m = {}", PARAMS.m);
}
