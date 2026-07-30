//! Cycle-accurate benchmarking of McEliece operations using RDTSC.
//!
//! # Usage
//!
//! Pin to a single core for consistent results:
//!
//! ```sh
//! taskset -c 0 cargo run --release --features mceliece460896 --example bench_cycles
//! ```
//!
//! # About RDTSC
//!
//! Reads the Time-Stamp Counter via the `rdtsc` instruction. On modern x86_64
//! CPUs this counts fixed-frequency reference cycles (invariant TSC), so it is
//! stable across frequency scaling. We use `lfence` before and after to prevent
//! out-of-order execution from skewing measurements.

use std::arch::x86_64::{_mm_lfence, _rdtsc};
use std::hint::black_box;

fn rdtsc() -> u64 {
    unsafe {
        _mm_lfence();
        let t = _rdtsc();
        _mm_lfence();
        t
    }
}

fn measure_overhead(samples: usize) -> Vec<u64> {
    let mut results = Vec::with_capacity(samples);
    for _ in 0..samples {
        let start = rdtsc();
        black_box(());
        let end = rdtsc();
        results.push(end.saturating_sub(start));
    }
    results
}

fn min(data: &[u64]) -> u64 {
    *data.iter().min().unwrap_or(&0)
}

fn median(data: &mut [u64]) -> u64 {
    data.sort_unstable();
    let mid = data.len() / 2;
    if data.len() % 2 == 0 {
        (data[mid - 1] + data[mid]) / 2
    } else {
        data[mid]
    }
}

fn mean(data: &[u64]) -> f64 {
    let sum: u64 = data.iter().sum();
    sum as f64 / data.len() as f64
}

/// Benchmark one function, returning raw cycle counts with overhead subtracted.
fn bench<F: Fn()>(f: F, samples: usize, warmup: usize, overhead_median: u64) -> Vec<u64> {
    for _ in 0..warmup {
        black_box(f());
    }

    let mut results = Vec::with_capacity(samples);
    for _ in 0..samples {
        let start = rdtsc();
        black_box(f());
        let end = rdtsc();
        results.push(end.saturating_sub(start).saturating_sub(overhead_median));
    }
    results
}

fn print_stats(label: &str, data: &[u64]) {
    println!("  {label}:");
    println!(
        "    min:    {:>10.0} cycles  ({:.2} Mcycles)",
        min(data),
        min(data) as f64 / 1_000_000.0
    );
    println!(
        "    median: {:>10.0} cycles  ({:.2} Mcycles)",
        median(&mut data.to_vec()),
        median(&mut data.to_vec()) as f64 / 1_000_000.0
    );
    println!(
        "    mean:   {:>10.0} cycles  ({:.2} Mcycles)",
        mean(data) as u64,
        mean(data) / 1_000_000.0
    );
}

fn main() {
    // Keygen is very slow (~10s/iteration), so use few samples.
    // Encaps/decaps are fast enough for more.
    let kg_samples = 5;
    let fast_samples = 30;
    let warmup = 3;

    // Measure RDTSC overhead
    let overhead = measure_overhead(fast_samples);
    let overhead_median = median(&mut overhead.clone());
    println!(
        "RDTSC overhead: min={} median={} mean={:.0} cycles\n",
        min(&overhead),
        overhead_median,
        mean(&overhead)
    );

    // ── Key Generation ──────────────────────────────────────────
    println!("--- Key Generation ({} samples) ---", kg_samples);
    let kg = bench(
        || {
            black_box(mceliece_rs::mceliece::keygen());
        },
        kg_samples,
        warmup,
        overhead_median,
    );
    print_stats("KeyGen", &kg);

    // Generate a key once for the remaining benchmarks
    let (pk, sk) = mceliece_rs::mceliece::keygen();

    // ── Encapsulation ───────────────────────────────────────────
    println!("\n--- Encapsulation ({} samples) ---", fast_samples);
    let enc = bench(
        || {
            black_box(mceliece_rs::mceliece::encapsulate(&pk));
        },
        fast_samples,
        warmup,
        overhead_median,
    );
    print_stats("Encaps", &enc);

    // ── Decapsulation (valid ciphertext) ────────────────────────
    let (ct, _) = mceliece_rs::mceliece::encapsulate(&pk);

    println!(
        "\n--- Decapsulation — valid CT ({} samples) ---",
        fast_samples
    );
    let dec = bench(
        || {
            black_box(mceliece_rs::mceliece::decapsulate(&ct, &sk));
        },
        fast_samples,
        warmup,
        overhead_median,
    );
    print_stats("Decaps (valid)", &dec);

    // ── Decapsulation (tampered ciphertext) ─────────────────────
    let mut ct_tampered = ct.clone();
    if !ct_tampered.is_empty() {
        ct_tampered[0] ^= 0xFF;
    }

    println!(
        "\n--- Decapsulation — tampered CT ({} samples) ---",
        fast_samples
    );
    let dec_tampered = bench(
        || {
            black_box(mceliece_rs::mceliece::decapsulate(&ct_tampered, &sk));
        },
        fast_samples,
        warmup,
        overhead_median,
    );
    print_stats("Decaps (tampered)", &dec_tampered);

    // ── Summary ─────────────────────────────────────────────────
    println!("\n=== Summary ===");
    println!("Operation              median (cycles)  est. ms @ 2.5 GHz");
    println!("{:-<55}", "");
    println!(
        "KeyGen                 {:>12.0}       {:>8.1} ms",
        median(&mut kg.clone()) as f64 / 1_000_000.0,
        median(&mut kg.clone()) as f64 / 2_500_000.0
    );
    println!(
        "Encapsulate            {:>12.0}       {:>8.3} ms",
        median(&mut enc.clone()) as f64 / 1_000_000.0,
        median(&mut enc.clone()) as f64 / 2_500_000.0
    );
    println!(
        "Decapsulate (valid)    {:>12.0}       {:>8.3} ms",
        median(&mut dec.clone()) as f64 / 1_000_000.0,
        median(&mut dec.clone()) as f64 / 2_500_000.0
    );
    println!(
        "Decapsulate (tampered) {:>12.0}       {:>8.3} ms",
        median(&mut dec_tampered.clone()) as f64 / 1_000_000.0,
        median(&mut dec_tampered.clone()) as f64 / 2_500_000.0
    );
}
