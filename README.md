# Constant-Time Classic McEliece in Rust

Yet another implementation of the Classic McEliece Post-Quantum Cryptography (PQC) algorithm, written in Rust. 

This repository focuses on eliminating timing side-channel vulnerabilities (like Chosen-Ciphertext Attacks) while strictly adhering to the NIST PQC specifications and Known Answer Tests (KAT).

## Project Layout
 
This is a Cargo workspace with two crates:
 
* [`core/`](./core) — the McEliece implementation itself (package name `mceliece_rs`). This is what you want if you're using the library from Rust.
* [`python_bindings/`](./python_bindings) — PyO3 bindings exposing `core` to Python as the `mceliece_rs` module. See that directory for build/usage instructions.

## Features
* **Constant-Time Execution:** End-to-end constant-time decapsulation. Conditional branching and early returns have been replaced with bitwise masking (`subtle::conditional_select`) and dummy operations to prevent timing leaks.
* **Implicit Rejection:** Secure handling of tampered ciphertexts to thwart CCA (Chosen-Ciphertext Attacks).
* **NIST Specification Compliant:** Retains the original rejection sampling in the fixed-weight vector generation to ensure 100% compatibility with official KAT vectors, while utilizing cache-safe assignments.
* **Memory Safety:** Built with Rust, guaranteeing protection against buffer overflows and memory leaks—critical for embedded cryptographic deployments.
* **Python Bindings:** Optional PyO3-based bindings (`python_bindings/`) let you call keygen/encapsulate/decapsulate from Python.

## Security Proof: Statistical Timing Analysis
To mathematically prove the absence of timing leaks at the hardware level, this implementation has been tested using the **dudect** (Dude, is my code constant time?) methodology.

A Welch's t-test was executed on a VPS, processing 1000 decapsulation iterations across valid and tampered ciphertexts. 

**Results:**
* `n`: +0.001M
* `max t`: -2.86153
* `max tau`: -0.09104
* `(5/tau)^2`: 3016

With a `max t` value well within the safe interval of `[-4.5, +4.5]`, the hardware-level constant-time execution is statistically verified.


## Current Status: Security & Correctness Over Performance
The current implementation intentionally prioritizes **cryptographic correctness, memory safety, and side-channel security** over raw execution speed. 

Presently, a full decapsulation cycle takes significant time compared to the NIST AVX2 implementations. This performance gap is expected and stems from:
1. **Strict Constant-Time Overhead:** The removal of all early-returns forces the CPU to evaluate the full Chien search and Patterson algorithm matrices regardless of the error state, inherently maximizing the cycle count to prevent timing leaks.
2. **Naive GF Arithmetic:** The Galois Field multiplication (`GF::mul`) currently relies on bit-by-bit carry-less multiplication without hardware acceleration.
3. **Lack of SIMD/Vectorization:** Operations are executed sequentially rather than utilizing CPU-specific vector extensions (like `PCLMULQDQ` or `AVX512`).

The next focus area will be performance profiling and optimization:
* Replacing naive GF arithmetic with bitsliced or hardware-accelerated polynomial multiplications.
* Integrating SIMD parallelism for matrix-vector operations.
* Exploring `#![no_std]` compliant assembly blocks for critical inner loops.

## Python Bindings
 
The `python_bindings/` crate exposes `core` to Python via [PyO3](https://pyo3.rs) + [maturin](https://www.maturin.rs/):
 
```bash
pip install maturin
cd python_bindings
maturin develop --release
```
 
```python
import mceliece_rs
 
pk, sk = mceliece_rs.keygen()
ciphertext, session_key = mceliece_rs.encapsulate(pk)
recovered = mceliece_rs.decapsulate(ciphertext, sk)
assert session_key == recovered
```

## Future Plans & Research: Streaming Public Key Architecture
The primary bottleneck for deploying Classic McEliece on resource-constrained embedded systems (like Cortex-M4 or ESP32 with ~256KB RAM) is its massive Public Key size (>1 MB). 

The next phase of this project introduces a **Streaming Architecture**:
* **Zero-Allocation `no_std` Core:** Porting the mathematical core to `#![no_std]` for bare-metal execution.
* **I/O Streaming PK:** Instead of loading the entire Public Key into RAM, the matrix $T$ will be streamed row-by-row via external flash directly into the CPU registers during the $C = e_1 \oplus (T \cdot e_2)$ matrix-vector multiplication.
* **RAM Optimization:** This approach will reduce the RAM footprint from Megabytes to mere Kilobytes, making high-security McEliece practical for smart cards and IoT devices, without introducing any new side-channel vulnerabilities (as the streamed matrix is public data).
