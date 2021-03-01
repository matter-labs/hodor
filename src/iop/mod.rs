use ff::PrimeField;

pub mod trivial_coset_combiner;
pub mod blake2s_trivial_iop;

pub mod optimized_iop;

/*

This module contains an IOP abstraction that is implied to be instantiated by the Merkle tree, but 
in principle any secure vector accumulator may be used for the same purpose

At the moment some trait constraints are over-restrictive and imply use of the hash functions that 
are PRF and return a byte array as a digest. In principle for an instantiation of the IOP as a Merkle tree just
collision resistant hash function is enough for interactive case and a single invocation of PRF is enough for Fiat-Shamir
heuristics (taking a Merkle tree root as a challenge for some other function)

*/

pub trait CosetInformation: Sized + Clone + Copy {
    const COSET_SIZE: usize;
}

pub trait CosetCombiner<F: PrimeField> {
    const EXPECTED_DEGREE: usize;
    const COSET_SIZE: usize;
    // type CosetData: CosetInformation;
    
    fn get_for_natural_index(leafs: &[F], natural_index: usize) -> &F;
    fn get_for_tree_index(leafs: &[F], tree_index: usize) -> &F;
    fn tree_index_into_natural_index(tree_index: usize) -> usize;
    fn natural_index_into_tree_index(natural_index: usize) -> usize;
    fn get_coset_for_natural_index(natural_index: usize, domain_size: usize) -> Vec<usize>;
    fn get_coset_for_tree_index(tree_index: usize, domain_size: usize) -> Vec<usize>;
}

pub trait LeafEncoder<F: PrimeField> {
    type Output;

    fn encode_leaf(value: &F) -> Self::Output;
}

pub trait HashEncoder<F: PrimeField> {
    type Input;

    fn interpret_hash(value: &Self::Input) -> F;
}

pub trait HashFunctionOutput: AsRef<[u8]> + Clone + Eq + PartialEq {}

pub trait IopTreeHasher<F: PrimeField> {
    type HashOutput: HashFunctionOutput;
    type LeafEncoder: LeafEncoder<F> + HashEncoder<F>; 

    fn hash_leaf(value: &F) -> Self::HashOutput;
    fn hash_encoded_leaf(value: &<Self::LeafEncoder as LeafEncoder<F>>::Output) -> Self::HashOutput;
    fn hash_node(values: &[Self::HashOutput], level: usize) -> Self::HashOutput;
}

pub trait IopTree<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Hasher: IopTreeHasher<F>;
    fn create(leafs: &[F]) -> Self;
    fn size(&self) -> u64;
    fn get_root(&self) -> <Self::Hasher as IopTreeHasher<F>>::HashOutput;
    fn encode_root_into_challenge(root: & <Self::Hasher as IopTreeHasher<F>>::HashOutput) -> F;
    fn get_challenge_scalar_from_root(&self) -> F;
    fn verify(root: &<Self::Hasher as IopTreeHasher<F>>::HashOutput, leaf_value: &F, path: &[<Self::Hasher as IopTreeHasher<F>>::HashOutput], index: usize) -> bool;
    fn get_path(&self, index: usize, leafs_values: &[F]) -> Vec< <Self::Hasher as IopTreeHasher<F>>::HashOutput >;
}

pub trait IopQuery<F: PrimeField>: 'static + PartialEq + Eq + Clone {
    type Hasher: IopTreeHasher<F>;

    fn tree_index(&self) -> usize;
    fn natural_index(&self) -> usize;
    fn value(&self) -> F;
    fn path(&self) ->  &[<Self::Hasher as IopTreeHasher<F>>::HashOutput];
}

pub trait IOP<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Tree: IopTree<F, Combiner = Self::Combiner>;
    type Query: IopQuery<F, Hasher = <Self::Tree as IopTree<F> >::Hasher>;

    fn create(leafs: & [F]) -> Self;
    fn get_for_natural_index(leafs: &[F], natural_index: usize) -> &F;
    fn get_for_tree_index(leafs: &[F], tree_index: usize) -> &F;
    fn get_root(&self) -> < <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput;
    fn encode_root_into_challenge(root: & < <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput) -> F;
    fn get_challenge_scalar_from_root(&self) -> F;
    fn verify_query(query: &Self::Query, root: &< <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput) -> bool;
    fn query(&self, natural_index: usize, leafs: &[F]) -> Self::Query;
}