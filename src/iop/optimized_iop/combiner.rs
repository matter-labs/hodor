// use crate::iop::CosetCombiner as InnerCosetCombiner;
// use ff::PrimeField;
// pub trait CosetCombiner<F: PrimeField> {
//     type InnerCombiner: InnerCosetCombiner<F>;    

//     fn get_for_natural_index(leafs: &[Vec<F>], natural_index: usize) -> Vec<F>;
//     fn get_for_tree_index(leafs: &[Vec<F>], tree_index: usize) -> Vec<F>;
//     fn tree_index_into_natural_index(tree_index: usize) -> usize;
//     fn natural_index_into_tree_index(natural_index: usize) -> usize;
//     fn get_coset_for_natural_index(natural_index: usize, domain_size: usize) -> Vec<usize>;
//     fn get_coset_for_tree_index(tree_index: usize, domain_size: usize) -> Vec<usize>;
// }

// pub struct TrivialCombiner<F: PrimeField> {
//     _marker: std::marker::PhantomData<F>,
// }

// impl<'c, F: PrimeField> CosetCombiner<F> for TrivialCombiner<F> {
//     type InnerCombiner = crate::iop::trivial_coset_combiner::TrivialCombiner<F>;

//     #[inline(always)]
//     fn get_for_natural_index(ldes: &[Vec<F>], natural_index: usize) -> Vec<F> {
//         let mut values = vec![F::zero(); ldes.len()];

//         for (idx, lde) in ldes.iter().enumerate() {
//             let value = <Self::InnerCombiner as InnerCosetCombiner<F>>::get_for_natural_index(
//                 &lde,
//                 natural_index,
//             );
//             values[idx] = *value;
//         }

//         values
//     }

//     #[inline(always)]
//     fn get_for_tree_index(ldes: &[Vec<F>], tree_index: usize) -> Vec<F> {
//         let mut values = vec![F::zero(); ldes.len()];

//         for (idx, lde) in ldes.iter().enumerate() {
//             let value = <Self::InnerCombiner as InnerCosetCombiner<F>>::get_for_tree_index(
//                 &lde, tree_index,
//             );
//             values[idx] = *value;
//         }

//         values
//     }

//     fn get_coset_for_natural_index(natural_index: usize, domain_size: usize) -> Vec<usize> {
//         assert!(natural_index < domain_size);
//         let natural_pair_index = (natural_index + (domain_size / 2)) % domain_size;
//         let mut coset = vec![natural_index, natural_pair_index];
//         coset.sort();

//         coset
//     }

//     fn get_coset_for_tree_index(natural_index: usize, domain_size: usize) -> Vec<usize> {
//         Self::get_coset_for_natural_index(natural_index, domain_size)
//     }

//     #[inline(always)]
//     fn tree_index_into_natural_index(tree_index: usize) -> usize {
//         tree_index
//     }

//     #[inline(always)]
//     fn natural_index_into_tree_index(natural_index: usize) -> usize {
//         natural_index
//     }
// }
