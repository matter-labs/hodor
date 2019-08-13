use ff::PrimeField;

pub mod trivial_coset_combiner;
pub mod blake2s_trivial_iop;

pub trait CosetCombiner<F: PrimeField> {
    type Index;
    fn shuffle_for_iop(values: Vec<F>) -> (Vec<F>, Vec<Self::Index>);
}

pub trait LeafEncoder<F: PrimeField> {
    type Output;

    fn encode_leaf(value: &F) -> Self::Output;
}

pub trait HashEncoder<F: PrimeField> {
    type Input;

    fn interpret_hash(value: &Self::Input) -> F;
}

pub trait IopTreeHasher<F: PrimeField> {
    type HashOutput: AsRef<[u8]>;
    type LeafEncoder: LeafEncoder<F> + HashEncoder<F>; 

    fn hash_leaf(value: &F) -> Self::HashOutput;
    fn hash_node(values: &[Self::HashOutput], level: usize) -> Self::HashOutput;
}

pub trait IopTree<F: PrimeField> {
    type Hasher: IopTreeHasher<F>;
    fn create(leafs: &[F]) -> Self;
    fn get_root(&self) -> <Self::Hasher as IopTreeHasher<F>>::HashOutput;
    fn get_challenge_scalar_from_root(&self) -> F;
}

pub trait IOP<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Tree: IopTree<F>;

    fn create(leafs: &[F]) -> Self;
    fn get_root(&self) -> < <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput;
    fn get_challenge_scalar_from_root(&self) -> F;
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow+1)) <= num {
        pow += 1;
    }

    pow
}