use ff::PrimeField;


pub trait CosetCombiner<F: PrimeField> {
    fn shuffle_for_iop(values: Vec<F>) -> Vec<F>;
}

pub trait LeafEncoder<F: PrimeField> {
    type Output;

    fn encode_leaf(value: &F) -> Self::Output;
}

pub trait IopTreeHasher<F: PrimeField> {
    type HashOutput;
    type LeafsEncoder: LeafEncoder<F>; 

    fn hash_leaf(value: &Self::LeafsEncoder::Output) -> Self::HashOutput;
    fn hash_node(values: &[Self::HashOutput], level: usize) -> Self::HashOutput;
}

pub trait IopTree<F: PrimeField> {
    type Hasher: IopTreeHasher<F>;
    fn create(leafs: &[F]) -> Self;
    // fn create_proof
}

pub trait IOP<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Tree: IopTree<F>;
}