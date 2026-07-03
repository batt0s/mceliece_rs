use dudect_bencher::{BenchRng, Class, CtRunner, ctbench_main};
use mceliece_rs::mceliece::{PrivateKey, PublicKey, decapsulate, encapsulate, seeded_keygen};
use rand::RngCore;
use std::sync::OnceLock;

static KEYPAIR: OnceLock<(PublicKey, PrivateKey, Vec<u8>)> = OnceLock::new();

fn bench_decapsulate(runner: &mut CtRunner, mut rng: &mut BenchRng) {
    let (_pk, sk, valid_c) = KEYPAIR.get_or_init(|| {
        let seed = [42u8; 32];
        let (pk, sk) = seeded_keygen(seed);
        let (c, _k) = encapsulate(&pk);
        (pk, sk, c)
    });

    let mut inputs = Vec::new();
    let mut classes = Vec::new();

    for _ in 0..1_000 {
        let is_tampered = rng.next_u32() % 2 == 0;
        let class = if is_tampered {
            Class::Right
        } else {
            Class::Left
        };

        let mut test_c = valid_c.clone();
        if is_tampered {
            test_c[0] ^= 0xFF; // Implicit rejection
        }

        inputs.push(test_c);
        classes.push(class);
    }

    for (class, test_c) in classes.into_iter().zip(inputs.into_iter()) {
        runner.run_one(class, || decapsulate(&test_c, sk));
    }
}

ctbench_main!(bench_decapsulate);
