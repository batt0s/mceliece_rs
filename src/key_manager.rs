use crate::mceliece::{PrivateKey, PublicKey, seeded_keygen};
use argon2::Argon2;
use chacha20poly1305::{
    ChaCha20Poly1305, Key, Nonce,
    aead::{Aead, KeyInit},
};
use rand::{RngCore, rngs::OsRng};
use std::fs;
use std::path::Path;

pub fn save_keys(
    public_key: &PublicKey,
    private_key: &PrivateKey,
    pub_path: &Path,
    priv_path: &Path,
    password: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let pk_bytes = bincode::serialize(public_key)?;
    fs::write(pub_path, pk_bytes)?;

    let mut salt = [0u8; 16];
    OsRng.fill_bytes(&mut salt);

    let mut key = [0u8; 32];
    Argon2::default()
        .hash_password_into(password.as_bytes(), &salt, &mut key)
        .map_err(|e| e.to_string())?;

    let cipher = ChaCha20Poly1305::new(Key::from_slice(&key));

    let mut nonce_bytes = [0u8; 12];
    OsRng.fill_bytes(&mut nonce_bytes);
    let nonce = Nonce::from_slice(&nonce_bytes);

    let ciphertext = cipher
        .encrypt(nonce, private_key.delta.as_ref())
        .map_err(|e| e.to_string())?;

    let mut priv_file_data = Vec::new();
    priv_file_data.extend_from_slice(&salt);
    priv_file_data.extend_from_slice(&nonce_bytes);
    priv_file_data.extend_from_slice(&ciphertext);

    fs::write(priv_path, priv_file_data)?;

    Ok(())
}

pub fn load_keys(
    pub_path: &Path,
    priv_path: &Path,
    password: &str,
) -> Result<(PublicKey, PrivateKey), Box<dyn std::error::Error>> {
    let pk_bytes = fs::read(pub_path)?;
    let loaded_pk: PublicKey = bincode::deserialize(&pk_bytes)?;

    let priv_data = fs::read(priv_path)?;
    if priv_data.len() < 16 + 12 + 32 {
        return Err("Invalid private key file".into());
    }
    let salt = &priv_data[0..16];
    let nonce = Nonce::from_slice(&priv_data[16..28]);
    let ciphertext = &priv_data[28..];

    let mut key = [0u8; 32];
    Argon2::default()
        .hash_password_into(password.as_bytes(), salt, &mut key)
        .map_err(|e| e.to_string())?;

    let cipher = ChaCha20Poly1305::new(Key::from_slice(&key));
    let decrypted = cipher
        .decrypt(nonce, ciphertext)
        .map_err(|e| e.to_string())?;

    if decrypted.len() != 32 {
        return Err("Invalid private key file".into());
    }

    let mut delta = [0u8; 32];
    delta.copy_from_slice(&decrypted);

    let (generated_pk, sk) = seeded_keygen(delta);

    if loaded_pk.T != generated_pk.T {
        return Err("Public key mismatch".into());
    }

    Ok((generated_pk, sk))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mceliece::keygen;
    use std::env;
    use std::fs;

    #[test]
    fn test_save_and_load_keys_success() {
        let (pk, sk) = keygen();
        let password = "super_secure_password_123";

        let temp_dir = env::temp_dir();
        let pub_path = temp_dir.join("test_mceliece_success.pub");
        let priv_path = temp_dir.join("test_mceliece_success");

        let save_result = save_keys(&pk, &sk, &pub_path, &priv_path, password);
        assert!(save_result.is_ok(), "Keys saved to disk successfully!");

        let load_result = load_keys(&pub_path, &priv_path, password);
        assert!(load_result.is_ok(), "Keys loaded from disk successfully!");

        let (loaded_pk, loaded_sk) = load_result.unwrap();

        assert_eq!(pk.T, loaded_pk.T, "Public key does not match!");
        assert_eq!(
            sk.delta, loaded_sk.delta,
            "Private key delta does not match!"
        );
        assert_eq!(sk.s, loaded_sk.s, "Private key 's' vector does not match!");

        let _ = fs::remove_file(pub_path);
        let _ = fs::remove_file(priv_path);
    }

    #[test]
    fn test_load_keys_wrong_password() {
        let (pk, sk) = keygen();
        let correct_password = "dogru_parola";
        let wrong_password = "yanlis_parola";

        let temp_dir = env::temp_dir();
        let pub_path = temp_dir.join("test_mceliece_wrong_pass.pub");
        let priv_path = temp_dir.join("test_mceliece_wrong_pass");

        let save_result = save_keys(&pk, &sk, &pub_path, &priv_path, correct_password);
        assert!(save_result.is_ok(), "Keys saved to disk successfully!");

        let load_result = load_keys(&pub_path, &priv_path, wrong_password);
        assert!(
            load_result.is_err(),
            "Security vulnerability: Wrong password accepted without error!"
        );

        let _ = fs::remove_file(pub_path);
        let _ = fs::remove_file(priv_path);
    }

    #[test]
    fn test_load_keys_corrupted_file() {
        let (pk, sk) = keygen();
        let password = "my_password";

        let temp_dir = env::temp_dir();
        let pub_path = temp_dir.join("test_mceliece_corrupt.pub");
        let priv_path = temp_dir.join("test_mceliece_corrupt");

        let save_result = save_keys(&pk, &sk, &pub_path, &priv_path, password);
        assert!(save_result.is_ok(), "Keys saved to disk successfully!");

        let mut priv_data = fs::read(&priv_path).unwrap();
        let last_idx = priv_data.len() - 1;
        priv_data[last_idx] ^= 0xFF; // Bitleri tersine çevir
        fs::write(&priv_path, priv_data).unwrap();

        let load_result = load_keys(&pub_path, &priv_path, password);
        assert!(
            load_result.is_err(),
            "Security vulnerability: Corrupted/played back file accepted without error!"
        );

        let _ = fs::remove_file(pub_path);
        let _ = fs::remove_file(priv_path);
    }
}
