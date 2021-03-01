use std::convert::TryInto;

use ff::PrimeField;

use crate::iop::{trivial_coset_combiner::TrivialCombiner, CosetCombiner, IopTreeHasher};

use super::{
    query::{IopQuery, TrivialBlake2sIopQuery},
    tree::{Blake2sIopTree, IopTree},
    tree_hasher::TreeHasher,
};

pub trait IOP<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Tree: IopTree<F, Combiner = Self::Combiner>;
    type Query: IopQuery<F>;

    fn create(ldes: &[F], stride: usize) -> Self;
    fn get_for_natural_index(leafs: &[F], natural_index: usize) -> &F;
    fn get_for_tree_index(leafs: &[F], tree_index: usize) -> &F;
    fn get_root(&self) -> <<<Self::Tree as IopTree<F> >::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput;
    fn encode_root_into_challenge(
        root: & <<<Self::Tree as IopTree<F> >::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
    ) -> F;
    fn get_challenge_scalar_from_root(&self) -> F;
    fn verify_query(
        query: &Self::Query,
        root: &<<<Self::Tree as IopTree<F> >::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
    ) -> bool;
    fn query(&self, natural_index: usize, leafs: &[F], stride: usize) -> Self::Query;
}

pub struct TrivialBlake2sIOP<F: PrimeField> {
    tree: Blake2sIopTree<F>,
}

impl<F: PrimeField> IOP<F> for TrivialBlake2sIOP<F> {
    type Combiner = TrivialCombiner<F>;

    type Tree = Blake2sIopTree<F>;

    type Query = TrivialBlake2sIopQuery<F>;

    fn create(ldes: &[F], stride: usize) -> Self {
        Self {
            tree: <Self::Tree as IopTree<F>>::create(ldes, stride),
        }
    }

    fn get_for_natural_index(leafs: &[F], natural_index: usize) -> &F {
        <Self::Combiner as CosetCombiner<F>>::get_for_natural_index(leafs, natural_index)
    }

    fn get_for_tree_index(leafs: &[F], tree_index: usize) -> &F {
        <Self::Combiner as CosetCombiner<F>>::get_for_tree_index(leafs, tree_index)
    }

    fn get_root(&self) -> <<<Self::Tree as IopTree<F> >::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput{
        self.tree.get_root()
    }

    fn encode_root_into_challenge(
        root: &<<<Self::Tree as IopTree<F> >::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
    ) -> F {
        <Self::Tree as IopTree<F>>::encode_root_into_challenge(root)
    }

    fn get_challenge_scalar_from_root(&self) -> F {
        self.tree.get_challenge_scalar_from_root()
    }

    fn verify_query(
        query: &Self::Query,
        root: &<<<Self::Tree as IopTree<F> >::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
    ) -> bool {
        Self::Tree::verify(root, &query.value(), &query.path(), query.tree_index())
    }

    fn query(&self, natural_index: usize, ldes: &[F], stride: usize) -> Self::Query {
        assert!(natural_index < self.tree.size() as usize);
        
        let length_of_single_lde = stride;        
        assert!(natural_index < length_of_single_lde);

        let tree_index =
            <Self::Combiner as CosetCombiner<F>>::natural_index_into_tree_index(natural_index);
        
        let path = self.tree.get_path(tree_index, ldes, stride);
        
        let number_of_ldes = ldes.len() / length_of_single_lde;
        let mut values = vec![F::zero(); number_of_ldes];
        let ldes_splitted: Vec<&[F]> = ldes.chunks(length_of_single_lde).collect();

        for (idx, lde) in ldes_splitted.iter().enumerate() {
            let value = lde[natural_index];
            values[idx] = value;
        }

        TrivialBlake2sIopQuery::<F> {
            index: natural_index,
            values: values,
            path: path,
        }
    }
}
