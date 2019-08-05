use tiny_keccak::Keccak;
use blake2s_simd;
use crypto::sha2::Sha256;
use crypto::digest::Digest;

pub mod merkle_tree;
pub mod poly;

pub trait Hasher: Sized + Clone {
    fn new(personalization: &[u8]) -> Self;
    fn update(&mut self, data: &[u8]);
    fn finalize(&mut self) -> Vec<u8>;
}

// #[derive(Clone)]
// pub struct BlakeHasher {
//     h: Blake2s
// }

// impl Hasher for BlakeHasher {
//     fn new(personalization: &[u8]) -> Self {
//         let h = Blake2s::with_params(32, &[], &[], personalization);

//         Self {
//             h: h
//         }
//     }

//     fn update(&mut self, data: &[u8]) {
//         self.h.update(data);
//     }

//     fn finalize(&mut self) -> Vec<u8> {
//         let new_h = Blake2s::with_params(32, &[], &[], &[]);
//         let h = std::mem::replace(&mut self.h, new_h);

//         let result = h.finalize();

//         result.as_ref().to_vec().clone()
//     }
// }

#[derive(Clone)]
pub struct Keccak256Hasher {
    h: Keccak
}

impl Hasher for Keccak256Hasher {
    fn new(personalization: &[u8]) -> Self {
        let mut h = Keccak::new_keccak256();
        h.update(personalization);

        Self {
            h: h
        }
    }

    fn update(&mut self, data: &[u8]) {
        self.h.update(data);
    }

    fn finalize(&mut self) -> Vec<u8> {
        let new_h = Keccak::new_keccak256();
        let h = std::mem::replace(&mut self.h, new_h);

        let mut res: [u8; 32] = [0; 32];
        h.finalize(&mut res);

        res[..].to_vec()
    }
}

#[derive(Clone)]
pub struct Sha256Hasher {
    h: Sha256
}

impl Hasher for Sha256Hasher {
    fn new(personalization: &[u8]) -> Self {
        let mut h = Sha256::new();
        h.input(personalization);

        Self {
            h: h
        }
    }

    fn update(&mut self, data: &[u8]) {
        self.h.input(data);
    }

    fn finalize(&mut self) -> Vec<u8> {
        let new_h = Sha256::new();
        let mut h = std::mem::replace(&mut self.h, new_h);

        let mut res: [u8; 32] = [0; 32];
        h.result(&mut res);

        res[..].to_vec()
    }
}

pub fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow+1)) <= num {
        pow += 1;
    }

    pow
}