use crate::utils::{Hasher};

pub trait TreeHasher: Hasher{
    fn hash_leaf(value: &[u8]) -> Vec<u8> {
        let mut hasher = Self::new(b"leaf");
        hasher.update(&value);
        let hash = hasher.finalize();

        hash
    }

    fn compress(lhs: &[u8], rhs: &[u8], _i: usize) -> Vec<u8> {
        let mut hasher = Self::new(&[]);
        hasher.update(&lhs[..]);
        hasher.update(&rhs[..]);
        let hash = hasher.finalize();

        hash

    }
}

pub mod sequential_smt;