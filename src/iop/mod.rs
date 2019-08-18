use ff::PrimeField;

pub mod trivial_coset_combiner;
pub mod blake2s_trivial_iop;

pub trait CosetInformation: Sized + Clone + Copy {
    const COSET_SIZE: usize;
}

pub trait CosetCombiner<F: PrimeField> {
    const EXPECTED_DEGREE: usize;
    const COSET_SIZE: usize;
    // type CosetData: CosetInformation;
    
    // fn new<'l>(leafs: &'l [F]) -> Self where 'l: 'c;

    // fn get(&self, natural_index: usize) -> &'c F;
    fn get_leaf(leafs: &[F], tree_index: usize) -> &F;
    fn get_coset_for_index(natural_index: usize, domain_size: usize) -> Vec<usize>;
    // fn get_coset_for_index_f(natural_index: usize, domain_size: usize) -> [usize; Self::COSET_SIZE];
    // fn shuffle_for_iop(values: Vec<F>) -> (Vec<F>, Vec<Self::Index>);
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

pub trait IopQuery<F: PrimeField>: 'static {
    type Hasher: IopTreeHasher<F>;

    fn index(&self) -> usize;
    fn value(&self) -> F;
    fn path(&self) ->  &[<Self::Hasher as IopTreeHasher<F>>::HashOutput];
}

pub trait IOP<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Tree: IopTree<F, Combiner = Self::Combiner>;
    type Query: IopQuery<F, Hasher = <Self::Tree as IopTree<F> >::Hasher>;

    fn create(leafs: & [F]) -> Self;
    fn get_combined(leafs: &[F], tree_index: usize) -> &F;
    fn get_root(&self) -> < <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput;
    fn encode_root_into_challenge(root: & < <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput) -> F;
    fn get_challenge_scalar_from_root(&self) -> F;
    fn verify_query(query: &Self::Query, root: &< <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput) -> bool;
    fn query(&self, index: usize, leafs: &[F]) -> Self::Query;
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow+1)) <= num {
        pow += 1;
    }

    pow
}