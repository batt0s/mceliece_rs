// You can get the KAT vectors from classic.mceliece.org: https://classic.mceliece.org/nist/mceliece-kat-20221023.tar.gz
use hex;
use mceliece_rs::mceliece::{encapsulate_with_rng, seeded_keygen};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::string::String;

mod common;
use common::NistAesCtrDrbg;

struct KatVector {
    seed: [u8; 48],
    pk: String,
    sk: String,
    ct: String,
    ss: String,
}

fn parse_kat(path: &str) -> Vec<KatVector> {
    let file = File::open(path).expect("Could not open KAT file!");
    let reader = BufReader::new(file);

    let mut vectors = Vec::new();
    let (mut seed, mut pk, mut sk, mut ct) =
        ([0u8; 48], String::new(), String::new(), String::new());

    for line_result in reader.lines() {
        let line = line_result.unwrap();
        if let Some(v) = line.strip_prefix("seed = ") {
            seed = hex::decode(v)
                .unwrap()
                .try_into()
                .expect("seed must be 48 bytes");
        } else if let Some(v) = line.strip_prefix("pk = ") {
            pk = v.to_string();
        } else if let Some(v) = line.strip_prefix("sk = ") {
            sk = v.to_string();
        } else if let Some(v) = line.strip_prefix("ct = ") {
            ct = v.to_string();
        } else if let Some(v) = line.strip_prefix("ss = ") {
            vectors.push(KatVector {
                seed,
                pk: pk.clone(),
                sk: sk.clone(),
                ct: ct.clone(),
                ss: v.to_string(),
            });
        }
    }
    vectors
}

fn run_kat_test(path: &str) {
    let vectors = parse_kat(path);
    assert!(!vectors.is_empty(), "No KAT vectors parsed");

    for (i, v) in vectors.iter().enumerate() {
        let mut drbg = NistAesCtrDrbg::new(&v.seed);
        let mut delta = [0u8; 32];
        drbg.randombytes(&mut delta);

        let (pk, sk) = seeded_keygen(delta);

        assert_eq!(
            hex::encode(pk.to_bytes()).to_uppercase(),
            v.pk.to_uppercase(),
            "Public Keys Does Not Match"
        );
        assert_eq!(
            hex::encode(sk.to_bytes()).to_uppercase(),
            v.sk.to_uppercase(),
            "Secret Key Does Not Match"
        );

        let (ct, ss) = encapsulate_with_rng(&pk, &mut drbg);

        assert_eq!(
            hex::encode(&ct).to_uppercase(),
            v.ct.to_uppercase(),
            "Ciphertext Does Not Match"
        );
        assert_eq!(
            hex::encode(&ss).to_uppercase(),
            v.ss.to_uppercase(),
            "Session Key Does Not Match"
        );
    }
}

#[cfg(feature = "mceliece348864")]
#[test]
fn test_kat_mceliece348864() {
    run_kat_test("tests/data/mceliece348864/kat_kem.rsp");
}

#[cfg(feature = "mceliece348864f")]
#[test]
fn test_kat_mceliece348864f() {
    run_kat_test("tests/data/mceliece348864f/kat_kem.rsp");
}

#[cfg(feature = "mceliece460896")]
#[test]
fn test_kat_mceliece460896() {
    run_kat_test("tests/data/mceliece460896/kat_kem.rsp");
}

#[cfg(feature = "mceliece460896f")]
#[test]
fn test_kat_mceliece460896f() {
    run_kat_test("tests/data/mceliece460896f/kat_kem.rsp");
}

#[cfg(feature = "mceliece6688128")]
#[test]
fn test_kat_mceliece6688128() {
    run_kat_test("tests/data/mceliece6688128/kat_kem.rsp");
}

#[cfg(feature = "mceliece6688128f")]
#[test]
fn test_kat_mceliece6688128f() {
    run_kat_test("tests/data/mceliece6688128f/kat_kem.rsp");
}

#[cfg(feature = "mceliece6960119")]
#[test]
fn test_kat_mceliece6960119() {
    run_kat_test("tests/data/mceliece6960119/kat_kem.rsp");
}

#[cfg(feature = "mceliece6960119f")]
#[test]
fn test_kat_mceliece6960119f() {
    run_kat_test("tests/data/mceliece6960119f/kat_kem.rsp");
}

#[cfg(feature = "mceliece8192128")]
#[test]
fn test_kat_mceliece8192128() {
    run_kat_test("tests/data/mceliece8192128/kat_kem.rsp");
}

#[cfg(feature = "mceliece8192128f")]
#[test]
fn test_kat_mceliece8192128f() {
    run_kat_test("tests/data/mceliece8192128f/kat_kem.rsp");
}
