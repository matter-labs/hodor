use ff::{PrimeField, PrimeFieldRepr};
use blake2s_simd::{Params, State};
use crate::fft::multicore::Worker;
use super::*;
use super::trivial_coset_combiner::*;

lazy_static! {
    static ref BASE_BLAKE2S_PARAMS: State = {
        Params::new()
            .hash_length(32)
            .key(b"Squeamish Ossifrage")
            .personal(b"Shaftoe")
            .to_state()
    };
}

pub struct Blake2sLeafEncoder<F: PrimeField> {
    _marker: std::marker::PhantomData<F>
}

impl<F: PrimeField> Blake2sLeafEncoder<F>  {
    const SHAVE_BITS: u32 = 256 - F::CAPACITY;
    pub fn new() -> Self {
        assert!(F::NUM_BITS < 256);

        Self {
            _marker: std::marker::PhantomData
        }
    }
}

impl<F: PrimeField> LeafEncoder<F> for Blake2sLeafEncoder<F>{
    type Output = [u8; 32];

    fn encode_leaf(value: &F) -> Self::Output {
        let raw_repr = value.into_raw_repr();
        let mut output = [0u8; 32];
        raw_repr.write_le(&mut output[..]).expect("will write");

        output
    }
}

impl<F: PrimeField> HashEncoder<F> for Blake2sLeafEncoder<F>{
    type Input = [u8; 32];

    fn interpret_hash(value: &Self::Input) -> F {
        let value = *value;
        let mut repr = F::Repr::default();
        let shaving_mask: u64 = 0xffffffffffffffff >> (Self::SHAVE_BITS % 64);
        repr.read_be(&value[..]).expect("will read");
        // repr.read_le(&value[..]).expect("will read");
        let last_limb_idx = repr.as_ref().len() - 1;
        repr.as_mut()[last_limb_idx] &= shaving_mask;
        // let value = F::from_raw_repr(repr).expect("in a field");
        let value = F::from_repr(repr).expect("in a field");

        value
    }
}

pub struct Blake2sTreeHasher<F: PrimeField> {
    encoder: Blake2sLeafEncoder<F>
}

impl<F: PrimeField> Blake2sTreeHasher<F> {
    pub fn new() -> Self {
        Self {
            encoder: Blake2sLeafEncoder::new()
        }
    }
}

impl<F: PrimeField> IopTreeHasher<F> for Blake2sTreeHasher<F> {
    type HashOutput = [u8; 32];
    type LeafEncoder = Blake2sLeafEncoder<F>;

    fn hash_leaf(value: &F) -> Self::HashOutput {
        let value = <Self::LeafEncoder as LeafEncoder<F> >::encode_leaf(value);
        Self::hash_encoded_leaf(&value)
    }

    fn hash_encoded_leaf(value: &<Self::LeafEncoder as LeafEncoder<F>>::Output) -> Self::HashOutput {
        let mut state = (*BASE_BLAKE2S_PARAMS).clone();
        state.update(value);
        let output = state.finalize();

        *output.as_array()
    }

    fn hash_node(values: &[Self::HashOutput], _level: usize) -> Self::HashOutput {
        debug_assert!(values.len() == 2);
        let mut state = (*BASE_BLAKE2S_PARAMS).clone();
        for value in values.iter() {
            state.update(value);
        }

        let output = state.finalize();

        *output.as_array()
    }
}

pub struct Blake2sIopTree<F: PrimeField> {
    size: u64,
    nodes: Vec< < Blake2sTreeHasher<F> as IopTreeHasher<F> >::HashOutput >,
    hasher: Blake2sTreeHasher<F>
}

impl<F: PrimeField> Blake2sIopTree<F> {
    pub fn new() -> Self {
        Self {
            size: 0u64,
            nodes: vec![],
            hasher: Blake2sTreeHasher::new()
        }
    }
}

impl<'a, F: PrimeField> IopTree<'a, F> for Blake2sIopTree<F> {
    type Combiner = TrivialCombiner<'a, F>;
    type Hasher = Blake2sTreeHasher<F>;

    fn size(&self) -> u64 {
        self.size
    }

    fn create<'l>(leafs: &'l [F]) -> Self where 'l: 'a {
        {
            let _ = *BASE_BLAKE2S_PARAMS;
        }

        let num_leafs = leafs.len();
        assert!(num_leafs == num_leafs.next_power_of_two());
        let num_nodes = num_leafs;

        let combiner = <Self::Combiner as CosetCombiner<'a, F>>::new(leafs);

        let size = num_leafs as u64;

        // TODO: May be leafs are redundant, decide later on better view on API

        // these are nodes and encoded (not hashed) leafs
        let mut nodes = vec![[0u8; 32]; num_nodes];

        let worker = Worker::new();

        let mut leaf_hashes = vec![[0u8; 32]; num_leafs];

        let combiner_ref = &combiner;

        {
            worker.scope(leafs.len(), |scope, chunk| {
                for (i, lh) in leaf_hashes.chunks_mut(chunk)
                                .enumerate() {
                    scope.spawn(move |_| {
                        let base_idx = i*chunk;
                        for (j, lh) in lh.iter_mut().enumerate() {
                            let idx = base_idx + j;
                            *lh = < Self::Hasher as IopTreeHasher<F> >::hash_leaf(combiner_ref.get(idx));
                        }
                    });
                }
            });
        }

        // leafs are now encoded and hashed, so let's make a tree

        let num_levels = log2_floor(num_leafs) as usize;
        let mut nodes_for_hashing = &mut nodes[..];

        // separately hash last level, which hashes leaf hashes into first nodes
        {
            let level = num_levels-1;
            let inputs = &mut leaf_hashes[..];
            let (_, outputs) = nodes_for_hashing.split_at_mut(nodes_for_hashing.len()/2);
            assert!(outputs.len() * 2 == inputs.len());
            assert!(outputs.len().is_power_of_two());

            worker.scope(outputs.len(), |scope, chunk| {
                for (o, i) in outputs.chunks_mut(chunk)
                                .zip(inputs.chunks(chunk*2)) {
                    scope.spawn(move |_| {
                        for (o, i) in o.iter_mut().zip(i.chunks(2)) {
                            // let input: Vec<_> = i.collect();
                            *o = < Self::Hasher as IopTreeHasher<F> >::hash_node(i, level);
                        }
                    });
                }
            });
        }

        for level in (0..(num_levels-1)).rev() {
            // do the trick - split
            let (next_levels, inputs) = nodes_for_hashing.split_at_mut(nodes_for_hashing.len()/2);
            let (_, outputs) = next_levels.split_at_mut(next_levels.len() / 2);
            assert!(outputs.len() * 2 == inputs.len());
            assert!(outputs.len().is_power_of_two());

            worker.scope(outputs.len(), |scope, chunk| {
                for (o, i) in outputs.chunks_mut(chunk)
                                .zip(inputs.chunks(chunk*2)) {
                    scope.spawn(move |_| {
                        for (o, i) in o.iter_mut().zip(i.chunks(2)) {
                            *o = < Self::Hasher as IopTreeHasher<F> >::hash_node(i, level);
                        }
                    });
                }
            });

            nodes_for_hashing = next_levels;
        }

        Self {
            size: size,
            nodes: nodes,
            hasher: Blake2sTreeHasher::new()
        }
    }

    fn get_root(&self) -> <Self::Hasher as IopTreeHasher<F>>::HashOutput {
        // 1 here is ok, cause we have nodes as power of two, but not 2^n - 1
        self.nodes[1]
    }

    fn encode_root_into_challenge(root: & <Self::Hasher as IopTreeHasher<F>>::HashOutput) -> F {
        < <Self::Hasher as IopTreeHasher<F> >::LeafEncoder as HashEncoder<F> >::interpret_hash(&root)
    }

    fn get_challenge_scalar_from_root(&self) -> F {
        let root = self.get_root();

        Self::encode_root_into_challenge(&root)
    }

    fn verify(root: &<Self::Hasher as IopTreeHasher<F>>::HashOutput, leaf_value: &F, path: &[<Self::Hasher as IopTreeHasher<F>>::HashOutput], index: usize) -> bool {
        let mut hash = <Self::Hasher as IopTreeHasher<F>>::hash_leaf(&leaf_value);
        let mut idx = index;
        for el in path.iter() {
            if idx & 1usize == 0 {
                hash = <Self::Hasher as IopTreeHasher<F>>::hash_node(&[hash, el.clone()], 0);
            } else {
                hash = <Self::Hasher as IopTreeHasher<F>>::hash_node(&[el.clone(), hash], 0);
            }
            idx >>= 1;
        }

        &hash == root
    }

    fn get_path<'l>(&self, index: usize, leafs_values: &'l [F]) -> Vec< <Self::Hasher as IopTreeHasher<F>>::HashOutput > where 'l: 'a {
        assert!(self.size == self.nodes.len() as u64);
        let mut nodes = &self.nodes[..];

        let pair_index = index ^ 1usize;

        let mut path = vec![];

        let pair = &leafs_values[pair_index as usize];
        let encoded_pair_hash = < Self::Hasher as IopTreeHasher<F> >::hash_leaf(pair);
        path.push(encoded_pair_hash);

        let mut idx = index as usize;
        idx >>= 1;

        for _ in 0..log2_floor(nodes.len() / 2) {
            let half_len = nodes.len() / 2;
            let (next_level, this_level) = nodes.split_at(half_len);
            let pair_idx = idx ^ 1usize;
            let value = this_level[pair_idx];
            path.push(value);
            idx >>= 1;
            nodes = next_level;
        }

        path
    }
}

pub struct TrivialBlake2sIOP<F: PrimeField> {
    tree: Blake2sIopTree<F>,
}


impl<'i, F: PrimeField> IOP<'i, F> for TrivialBlake2sIOP<F> {
    type Combiner = TrivialCombiner<'i, F>;
    type Tree = Blake2sIopTree<F>;
    type Query = TrivialBlake2sIopQuery<F>;

    fn create<'l> (leafs: &'l [F]) -> Self{
        let tree = Self::Tree::create(leafs);

        Self {
            tree
        }
    }

    fn combine<'l>(leafs: &'l [F]) -> Self::Combiner where 'l: 'i {
        <Self::Combiner as CosetCombiner<'i, F>>::new(leafs)
    }

    fn get_root(&self) -> < <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput {
        self.tree.get_root()
    }

    fn encode_root_into_challenge(root: & < <Self::Tree as IopTree<'i, F> >::Hasher as IopTreeHasher<F>>::HashOutput) -> F {
        <Self::Tree as IopTree<'i, F> >::encode_root_into_challenge(root)
    }

    fn get_challenge_scalar_from_root(&self) -> F {
        self.tree.get_challenge_scalar_from_root()
    }

    fn verify_query(query: &Self::Query, root: &< <Self::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput) -> bool {
        Self::Tree::verify(root, &query.value(), &query.path(), query.index())
    }

    fn query(&self, index: usize, leafs: &[F]) -> Self::Query {
        assert!(index < self.tree.size() as usize);
        assert!(index < leafs.len());
        let value = leafs[index as usize];

        let path = self.tree.get_path(index, leafs);

        TrivialBlake2sIopQuery::<F> {
            index: index,
            value: value,
            path: path
        }
    }
}

#[derive(Clone, Debug)]
pub struct TrivialBlake2sIopQuery<F: PrimeField> {
    index: usize,
    value: F,
    path: Vec<[u8; 32]>,
}

impl<F: PrimeField> IopQuery<F> for TrivialBlake2sIopQuery<F> {
    type Hasher = Blake2sTreeHasher<F>;

    fn index(&self) -> usize {
        self.index
    }

    fn value(&self) -> F {
        self.value
    }

    fn path(&self) ->  &[<Self::Hasher as IopTreeHasher<F>>::HashOutput] {
        &self.path
    }
}

#[test]
fn make_small_tree() {
    use ff::Field;
    use hex::encode;
    use crate::experiments::Fr;
    let inputs = vec![Fr::one(); 16];

    let tree = Blake2sIopTree::create(&inputs);
    let root = tree.get_root();
    let challenge_scalar = tree.get_challenge_scalar_from_root();
    println!("Root = {}, scalar = {}", encode(&root), challenge_scalar);
}

#[test]
fn make_small_iop() {
    use crate::iop::IOP;
    use ff::Field;
    const SIZE: usize = 4;
    use crate::experiments::Fr;
    let mut inputs = vec![];
    let mut f = Fr::one();
    for _ in 0..SIZE {
        inputs.push(f);
        f.double();
    }

    let iop = TrivialBlake2sIOP::create(&inputs);
    let root = iop.get_root();
    for i in 0..SIZE {
        let query = iop.query(i, &inputs);
        let valid = TrivialBlake2sIOP::verify_query(&query, &root);
        assert!(valid, "invalid query for leaf {}", i);
    }
}