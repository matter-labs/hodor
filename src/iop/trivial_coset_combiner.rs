use ff::PrimeField;
use super::*;

pub struct TrivialCombiner;

impl<F: PrimeField> CosetCombiner<F> for TrivialCombiner {
    type Index = ();

    #[inline(always)] 
    fn shuffle_for_iop(values: Vec<F>) -> (Vec<F>, Vec<Self::Index>) {
        (values, vec![])
    }
}