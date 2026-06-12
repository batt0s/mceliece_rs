use mceliece::key_manager;
use mceliece::mceliece::keygen;
use std::path::Path;

fn main() {
    let pub_file = Path::new("pub.pem");
    let priv_file = Path::new("priv.pem");
    let passwd = "super_secret_password";

    println!("Generating new key pair");
    let (pk, sk) = keygen();

    println!("Saving key pair to disk");
    key_manager::save_keys(&pk, &sk, pub_file, priv_file, passwd).unwrap();

    println!("Key pair generated and saved successfully");

    println!("Loading key pair from disk");
    let (pk, sk) = key_manager::load_keys(pub_file, priv_file, passwd).unwrap();

    println!("Key pair loaded successfully");
}
