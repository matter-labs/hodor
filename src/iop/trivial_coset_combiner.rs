use ff::PrimeField;
use super::*;

pub struct TrivialCombiner<'c, F: PrimeField> {
    leafs: & 'c [F]
}

impl<'c, F: PrimeField> CosetCombiner<'c, F> for TrivialCombiner<'c, F> {
    #[inline(always)] 
    fn get(&self, natural_index: usize) -> &'c F {
        &self.leafs[natural_index]
    }

    fn new<'l>(leafs: &'l [F]) -> Self where 'l: 'c {
        Self {
            leafs: leafs
        }
    }

    
    // fn shuffle_for_iop(values: Vec<F>) -> (Vec<F>, Vec<Self::Index>) {
    //     (values, vec![])
    // }
}