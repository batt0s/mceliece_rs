use mceliece_rs::key_manager;
use mceliece_rs::mceliece as mc;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::path::Path;

fn to_py_err<E: std::fmt::Display>(err: E) -> PyErr {
    PyValueError::new_err(err.to_string())
}

#[pyclass(name = "PublicKey", module = "mceliece_py")]
pub struct PyPublicKey(pub(crate) mc::PublicKey);

#[pymethods]
impl PyPublicKey {
    /// Serialize the public key (bincode-encoded `T` matrix).
    fn to_bytes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let bytes = bincode::serialize(&self.0).map_err(to_py_err)?;
        Ok(PyBytes::new(py, &bytes))
    }

    /// Deserialize a public key previously produced by `to_bytes()`.
    #[staticmethod]
    fn from_bytes(data: &[u8]) -> PyResult<Self> {
        let pk: mc::PublicKey = bincode::deserialize(data).map_err(to_py_err)?;
        Ok(PyPublicKey(pk))
    }

    fn __repr__(&self) -> String {
        format!("PublicKey(T_len={} u64 words)", self.0.T.len())
    }
}

#[pyclass(name = "PrivateKey", module = "mceliece_py")]
pub struct PyPrivateKey(pub(crate) mc::PrivateKey);

#[pymethods]
impl PyPrivateKey {
    /// The 32-byte seed this private key (and its matching public key) were
    /// deterministically derived from.
    fn delta<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.0.delta)
    }

    fn __repr__(&self) -> String {
        "PrivateKey(<redacted>)".to_string()
    }
}

#[pyfunction]
fn keygen(py: Python<'_>) -> PyResult<(PyPublicKey, PyPrivateKey)> {
    let (pk, sk) = py.detach(mc::keygen);
    Ok((PyPublicKey(pk), PyPrivateKey(sk)))
}

#[pyfunction]
fn seeded_keygen(py: Python<'_>, seed: &[u8]) -> PyResult<(PyPublicKey, PyPrivateKey)> {
    let seed_arr: [u8; 32] = seed.try_into().map_err(|_| {
        PyValueError::new_err(format!("seed must be exactly 32 bytes, got {}", seed.len()))
    })?;
    let (pk, sk) = py.detach(move || mc::seeded_keygen(seed_arr));
    Ok((PyPublicKey(pk), PyPrivateKey(sk)))
}

#[pyfunction]
fn encapsulate<'py>(
    py: Python<'py>,
    pk: &PyPublicKey,
) -> PyResult<(Bound<'py, PyBytes>, Bound<'py, PyBytes>)> {
    let (ciphertext, session_key) = py.detach(|| mc::encapsulate(&pk.0));
    Ok((
        PyBytes::new(py, &ciphertext),
        PyBytes::new(py, &session_key),
    ))
}

#[pyfunction]
fn decapsulate<'py>(py: Python<'py>, ciphertext: &[u8], sk: &PyPrivateKey) -> Bound<'py, PyBytes> {
    let ciphertext = ciphertext.to_vec();
    let key = py.detach(|| mc::decapsulate(&ciphertext, &sk.0));
    PyBytes::new(py, &key)
}

#[pyfunction]
fn save_keys(
    pk: &PyPublicKey,
    sk: &PyPrivateKey,
    pub_path: &str,
    priv_path: &str,
    password: &str,
) -> PyResult<()> {
    key_manager::save_keys(
        &pk.0,
        &sk.0,
        Path::new(pub_path),
        Path::new(priv_path),
        password,
    )
    .map_err(|e| PyIOError::new_err(e.to_string()))
}

#[pyfunction]
fn load_keys(
    pub_path: &str,
    priv_path: &str,
    password: &str,
) -> PyResult<(PyPublicKey, PyPrivateKey)> {
    let (pk, sk) = key_manager::load_keys(Path::new(pub_path), Path::new(priv_path), password)
        .map_err(|e| PyIOError::new_err(e.to_string()))?;
    Ok((PyPublicKey(pk), PyPrivateKey(sk)))
}

#[pymodule]
fn mceliece_py(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPublicKey>()?;
    m.add_class::<PyPrivateKey>()?;
    m.add_function(wrap_pyfunction!(keygen, m)?)?;
    m.add_function(wrap_pyfunction!(seeded_keygen, m)?)?;
    m.add_function(wrap_pyfunction!(encapsulate, m)?)?;
    m.add_function(wrap_pyfunction!(decapsulate, m)?)?;
    m.add_function(wrap_pyfunction!(save_keys, m)?)?;
    m.add_function(wrap_pyfunction!(load_keys, m)?)?;
    Ok(())
}
