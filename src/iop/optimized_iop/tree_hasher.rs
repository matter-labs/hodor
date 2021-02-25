use ff::PrimeField;

use crate::iop::{IopTreeHasher as InnerIopTreeHasher, LeafEncoder, blake2s_trivial_iop::Blake2sTreeHasher as InnerBlake2sTreeHasher};

pub trait TreeHasher<F: PrimeField> {
    type InnerIopTreeHasher: InnerIopTreeHasher<F>;

    fn hash_leaf(values: &[F]) -> <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput;
    // TODO: LeafEncoder + HashEncoder?
    fn hash_encoded_leaf(
        value: &<<Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::LeafEncoder as LeafEncoder<F>>::Output,
    ) -> <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput;
    fn hash_node(
        values: &[<Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput],
        level: usize,
    ) -> <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput;
}


pub struct Blake2sTreeHasher<F: PrimeField> {
    _m: std::marker::PhantomData<F>,
}

impl<F: PrimeField> Blake2sTreeHasher<F>{
    pub fn new() -> Self{
        Self{
            _m: std::marker::PhantomData,
        }
    }
}

impl<F: PrimeField> TreeHasher<F> for Blake2sTreeHasher<F> {
    type InnerIopTreeHasher = InnerBlake2sTreeHasher<F>;

    fn hash_leaf(values: &[F]) -> <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput {
        let mut buf = vec![[0u8; 32]; values.len()];
        for (idx, value) in values.iter().enumerate() {
            let encoded_value = <<Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::LeafEncoder as LeafEncoder<
                F,
            >>::encode_leaf(value);
            buf[idx].copy_from_slice(&encoded_value);
        }

        <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::hash_node(&buf, 0)
    }

    fn hash_encoded_leaf(
        value: &<<Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::LeafEncoder as LeafEncoder<F>>::Output,
    ) -> <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput {
        <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::hash_encoded_leaf(&value)
    }

    fn hash_node(
        values: &[<Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput],
        _level: usize,
    ) -> <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::HashOutput {
        <Self::InnerIopTreeHasher as InnerIopTreeHasher<F>>::hash_node(&values, _level)
    }
}
