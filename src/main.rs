use mceliece_rs::key_manager;
use mceliece_rs::mceliece::keygen;
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

#[cfg(test)]
mod ct_tests {
    use subtle::{Choice, ConditionallySelectable};

    // VURGUN YAPACAĞIMIZ NORMAL FONKSİYON (Zamanlaması tehlikeli)
    fn normal_secim(gizli_deger: u8, secenek_a: u8, secenek_b: u8) -> u8 {
        if gizli_deger == 42 {
            secenek_a
        } else {
            secenek_b
        }
    }

    // AYNI İŞİ YAPAN CONSTANT-TIME FONKSİYON (Güvenli)
    fn ct_secim(gizli_deger: u8, secenek_a: u8, secenek_b: u8) -> u8 {
        // 1. Karşılaştırma yapıyoruz ama "bool" yerine "Choice" dönüyor.
        // ct_eq (constant time equals) metodu subtle'dan gelir.
        // Eğer değerler eşitse Choice(1), değilse Choice(0) döner.
        let esit_mi: Choice = subtle::ConstantTimeEq::ct_eq(&gizli_deger, &42);

        // 2. if-else yerine seçimi yapıyoruz.
        // Eğer esit_mi == Choice(0) ise secenek_b'yi (ikinci parametreyi if'in else'i gibi düşün) döner.
        // Eğer esit_mi == Choice(1) ise secenek_a'yı döner.
        u8::conditional_select(&secenek_b, &secenek_a, esit_mi)
    }

    #[test]
    fn test_ikisi_de_ayni_calisir() {
        // Eşit olma durumu
        assert_eq!(normal_secim(42, 100, 200), 100);
        assert_eq!(ct_secim(42, 100, 200), 100);

        // Eşit olmama durumu
        assert_eq!(normal_secim(10, 100, 200), 200);
        assert_eq!(ct_secim(10, 100, 200), 200);
    }
}
