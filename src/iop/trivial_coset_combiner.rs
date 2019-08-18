use ff::PrimeField;
use super::*;

#[derive(Copy, Clone)]
pub struct CosetOfSizeTwo;

impl CosetInformation for CosetOfSizeTwo {
    const COSET_SIZE: usize = 2usize;
}

pub struct TrivialCombiner<'c, F: PrimeField> {
    leafs: & 'c [F]
}



impl<'c, F: PrimeField> CosetCombiner<'c, F> for TrivialCombiner<'c, F> {
    const EXPECTED_DEGREE: usize = 2usize;
    const COSET_SIZE: usize = 2usize;

    #[inline(always)] 
    fn get(&self, natural_index: usize) -> &'c F {
        &self.leafs[natural_index]
    }

    fn new<'l>(leafs: &'l [F]) -> Self where 'l: 'c {
        Self {
            leafs: leafs
        }
    }

    fn get_coset_for_index(natural_index: usize, domain_size: usize) -> Vec<usize> {
        assert!(natural_index < domain_size);
        let natural_pair_index = (natural_index + (domain_size / 2)) % domain_size;
        let mut coset = vec![natural_index, natural_pair_index];
        coset.sort();

        coset
    }
    
    // fn shuffle_for_iop(values: Vec<F>) -> (Vec<F>, Vec<Self::Index>) {
    //     (values, vec![])
    // }
}