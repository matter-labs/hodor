use ff::PrimeField;
use crate::SynthesisError;
use crate::domains::Domain;
use crate::polynomials::*;

use super::*;
use crate::iop::*;
use crate::iop::optimized_iop::IOP;

impl<'a, F: PrimeField, I: IOP<F>> FRIProofPrototype<F, I> {
    pub fn produce_proof(
        self,
        iop_values: &Polynomial<F, Values>,
        natural_first_element_index: usize,
    ) -> Result<FRIProof<F, I>, SynthesisError> {
        let mut domain_size = self.initial_degree_plus_one * self.lde_factor;
        let mut domain_idx = natural_first_element_index;

        let mut queries = vec![];
        let mut roots = vec![];

        for (iop, leaf_values) in Some(self.l0_commitment).iter().chain(&self.intermediate_commitments)
                                    .zip(Some(iop_values).into_iter().chain(&self.intermediate_values)) {
            
            let coset_values = <I::Combiner as CosetCombiner<F>>::get_coset_for_natural_index(domain_idx, domain_size);

            if coset_values.len() != <I::Combiner as CosetCombiner<F>>::COSET_SIZE {
                return Err(SynthesisError::InvalidValue(format!("invalid coset size, expected {}, got {}", <I::Combiner as CosetCombiner<F>>::COSET_SIZE, coset_values.len())));
            }

            for idx in coset_values.into_iter() {
                let query = iop.query(idx, leaf_values.as_ref(), leaf_values.as_ref().len());
                queries.push(query);
            }

            roots.push(iop.get_root());

            let (next_idx, next_size) = Domain::<F>::index_and_size_for_next_domain(domain_idx, domain_size);

            domain_idx = next_idx;
            domain_size = next_size;
        }

        let proof = FRIProof::<F, I> {
            queries,
            roots,
            final_coefficients: self.final_coefficients,
            initial_degree_plus_one: self.initial_degree_plus_one,
            output_coeffs_at_degree_plus_one: self.output_coeffs_at_degree_plus_one,
            lde_factor: self.lde_factor,
        };

        Ok(proof)
    }
}