use ff::PrimeField;

use crate::iop::{CosetCombiner, HashEncoder, IopTreeHasher, trivial_coset_combiner::TrivialCombiner};

use super::{
    tree_hasher::{Blake2sTreeHasher, TreeHasher},
};

use crate::utils::log2_floor;

use crate::fft::multicore::Worker;

use crate::iop::blake2s_trivial_iop::BASE_BLAKE2S_PARAMS;
pub trait IopTree<F: PrimeField> {
    type Combiner: CosetCombiner<F>;
    type Hasher: TreeHasher<F>;

    fn create(leafs: &[F], stride: usize) -> Self;
    fn size(&self) -> u64;
    fn get_root(
        &self,
    ) -> <<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput;
    fn encode_root_into_challenge(
        root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
    ) -> F;
    fn get_challenge_scalar_from_root(&self) -> F;
    fn verify(
        root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
        leaf_value: &[F],
        path: &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput],
        index: usize,
    ) -> bool;
    fn get_path(
        &self,
        index: usize,
        ldes: &[F],
        stride: usize,
    ) -> Vec<<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput>;
}

pub struct Blake2sIopTree<F: PrimeField> {
    size: u64,
    nodes: Vec<<<Blake2sTreeHasher<F> as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput>,
    hasher: Blake2sTreeHasher<F>,
}

impl<F: PrimeField> IopTree<F> for Blake2sIopTree<F> {    
    type Combiner = TrivialCombiner<F>;

    type Hasher = super::tree_hasher::Blake2sTreeHasher<F>;

    fn create(ldes: &[F], stride: usize) -> Self {
        
        {
            let _ = *BASE_BLAKE2S_PARAMS;
        }
        assert!(!ldes.is_empty());

        let num_ldes = ldes.len() / stride;

        let num_leafs = stride;

        let splitted_ldes : Vec<&[F]> = ldes.chunks_exact(stride).collect();
        for lde in splitted_ldes.iter() {
            assert_eq!(num_leafs, lde.len());
        }
        assert!(num_leafs == num_leafs.next_power_of_two());
        let num_nodes = num_leafs;

        let size = num_leafs as u64;

        // TODO: May be leafs are redundant, decide later on better view on API

        // these are nodes and encoded (not hashed) leafs
        let mut nodes = vec![[0u8; 32]; num_nodes];

        let worker = Worker::new();

        
        let mut leaf_hashes = vec![[0u8; 32]; num_leafs];
        
        let splitted_ldes_as_ref: &[&[F]] = splitted_ldes.as_ref();

        {
            worker.scope(num_leafs, |scope, chunk| {
                for (i, lh) in leaf_hashes.chunks_mut(chunk)
                                .enumerate() {
                    scope.spawn(move |_| {
                        let base_idx = i*chunk;
                        for (j, lh) in lh.iter_mut().enumerate() {
                            let idx = base_idx + j;
                            // collect each corresponding leafs
                            let mut values = vec![F::zero(); num_ldes];
                            for (lde_idx, lde) in splitted_ldes_as_ref.iter().enumerate(){
                                let value = <Self::Combiner as CosetCombiner<F> >::get_for_tree_index(&lde, idx);
                                values[lde_idx] = value.clone();
                            }

                            *lh = <Self::Hasher as TreeHasher<F>>::hash_leaf(&values);                            
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
            let level = num_levels - 1;
            let inputs = &mut leaf_hashes[..];
            let (_, outputs) = nodes_for_hashing.split_at_mut(nodes_for_hashing.len() / 2);
            assert!(outputs.len() * 2 == inputs.len());
            assert!(outputs.len().is_power_of_two());

            worker.scope(outputs.len(), |scope, chunk| {
                for (o, i) in outputs.chunks_mut(chunk).zip(inputs.chunks(chunk * 2)) {
                    scope.spawn(move |_| {
                        for (o, i) in o.iter_mut().zip(i.chunks(2)) {
                            // let input: Vec<_> = i.collect();
                            *o = <Self::Hasher as TreeHasher<F>>::hash_node(i, level);
                        }
                    });
                }
            });
        }


        for level in (0..(num_levels - 1)).rev() {
            // do the trick - split
            let (next_levels, inputs) = nodes_for_hashing.split_at_mut(nodes_for_hashing.len() / 2);
            let (_, outputs) = next_levels.split_at_mut(next_levels.len() / 2);
            assert!(outputs.len() * 2 == inputs.len());
            assert!(outputs.len().is_power_of_two());

            worker.scope(outputs.len(), |scope, chunk| {
                for (o, i) in outputs.chunks_mut(chunk).zip(inputs.chunks(chunk * 2)) {
                    scope.spawn(move |_| {
                        for (o, i) in o.iter_mut().zip(i.chunks(2)) {
                            *o = <Self::Hasher as TreeHasher<F>>::hash_node(i, level);
                        }
                    });
                }
            });

            nodes_for_hashing = next_levels;
        }


        Self {
            size: size,
            nodes: nodes,
            hasher: Blake2sTreeHasher::new(),
        }
    }

    fn size(&self) -> u64 {
        self.size
    }

    fn get_root(
        &self,
    ) -> <<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput {
        self.nodes[1]
    }

    fn encode_root_into_challenge(
        root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
    ) -> F {
        <<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::LeafEncoder as HashEncoder<F>>::interpret_hash(&root)
    }

    fn get_challenge_scalar_from_root(&self) -> F {
        Self::encode_root_into_challenge(&self.get_root())
    }

    fn verify(
        root: &<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput,
        leaf_values: &[F],
        path: &[<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput],
        tree_index: usize,
    ) -> bool {
        let mut hash = <Self::Hasher as TreeHasher<F>>::hash_leaf(&leaf_values);
        let mut idx = tree_index;
        for el in path.iter() {
            if idx & 1usize == 0 {
                hash = <Self::Hasher as TreeHasher<F>>::hash_node(&[hash, el.clone()], 0);
            } else {
                hash = <Self::Hasher as TreeHasher<F>>::hash_node(&[el.clone(), hash], 0);
            }
            idx >>= 1;
        }        
        &hash == root
    }

    // this function accepts tree index instead of natural index
    // thats why we need leaf values here
    fn get_path(
        &self,
        tree_index: usize,
        ldes: &[F],
        stride: usize,
    ) -> Vec<<<Self::Hasher as TreeHasher<F>>::InnerIopTreeHasher as IopTreeHasher<F>>::HashOutput>
    {
        assert!(self.size == self.nodes.len() as u64);
        let mut nodes = &self.nodes[..];

        let tree_pair_index = tree_index ^ 1usize;

        let mut path = vec![];

        let pair_natural_index =
            <Self::Combiner as CosetCombiner<F>>::tree_index_into_natural_index(tree_pair_index);        

        let number_of_ldes = ldes.len() / stride;

        let mut pair_values = vec![F::zero(); number_of_ldes];

        let ldes_splitted: Vec<&[F]> = ldes.chunks(stride).collect();
        for (idx, lde) in ldes_splitted.iter().enumerate(){
            pair_values[idx]   = lde[pair_natural_index];
        }

        // TODO 
        // for idx in 0..number_of_ldes
        //      pair_values[idx] = lde[idx*stride+pair_natural_index]

        let encoded_pair_hash = <Self::Hasher as TreeHasher<F>>::hash_leaf(&pair_values);
        path.push(encoded_pair_hash);

        let mut idx = tree_index;
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
