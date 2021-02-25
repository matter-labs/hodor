use ff::PrimeField;

use super::tree_hasher::{Blake2sTreeHasher, TreeHasher};
use crate::iop::IopTreeHasher;

pub trait IopQuery<F: PrimeField>: 'static + PartialEq + Eq + Clone {
    type Hasher: TreeHasher<F>;

    fn tree_index(&self) -> usize;
    fn natural_index(&self) -> usize;
    fn value(&self) -> &[F];
    fn path(
        &self,
    ) -> &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput];

    fn verify(
        root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
        leaf_values: &[F],
        path: &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput],
        tree_index: usize,
    );
}


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TrivialBlake2sIopQuery<F: PrimeField> {
    pub index: usize,
    pub values: Vec<F>,
    pub path: Vec<[u8; 32]>,
}

impl<F: PrimeField> IopQuery<F> for TrivialBlake2sIopQuery<F> {
    type Hasher = Blake2sTreeHasher<F>;

    fn natural_index(&self) -> usize {
        self.index
    }

    fn tree_index(&self) -> usize {
        self.index
    }

    fn value(&self) -> &[F] {
        &self.values[..] // TODO
    }

    fn path(&self) -> &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput] {
        &self.path
    }

    fn verify(
        root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
        leaf_values: &[F],
        path: &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput],
        tree_index: usize,
    ){
        todo!();
    }
}
// #[derive(Clone, Debug, PartialEq, Eq)]
// pub struct TrivialBlake2sIopQueryForSingle<F: PrimeField> {
//     pub index: usize,
//     pub value: F,
//     pub path: Vec<[u8; 32]>,
// }

// impl<F: PrimeField> IopQuery<F> for TrivialBlake2sIopQueryForSingle<F> {
//     type Hasher = Blake2sTreeHasher<F>;
//     type Value = F;

//     fn natural_index(&self) -> usize {
//         self.index
//     }

//     fn tree_index(&self) -> usize {
//         self.index
//     }

//     fn values(&self) -> Self::Value {
//         self.value.clone() // TODO
//     }

//     fn path(&self) -> &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput] {
//         &self.path
//     }


//     fn verify(
//         root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
//         leaf_values: &[F],
//         path: &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput],
//         tree_index: usize,
//     ){
//         todo!();
//     }
// }