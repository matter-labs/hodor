use ff::PrimeField;
use crate::SynthesisError;
use crate::domains::Domain;
use crate::polynomials::*;

use super::*;
use crate::iop::*;

impl<'a, F: PrimeField, I: IOP<F>> FRIProofPrototype<F, I> {
    // TODO: This is a draft and will transform into query checker after debugging
    pub fn produce_proof(
        self,
        iop_values: &Polynomial<F, Values>,
        natural_first_element_index: usize,
    ) -> Result<FRIProof<F, I>, SynthesisError> {
        let mut domain_size = self.initial_degree_plus_one * self.lde_factor;
        let mut next_domain_idx = natural_first_element_index;

        let mut queries = vec![];
        let mut roots = vec![];

        for (iop, leaf_values) in Some(self.l0_commitment).iter().chain(&self.intermediate_commitments)
                                    .zip(Some(iop_values).into_iter().chain(&self.intermediate_values)) {
            
            let coset_values = <I::Combiner as CosetCombiner<F>>::get_coset_for_natural_index(next_domain_idx, domain_size);

            if coset_values.len() != <I::Combiner as CosetCombiner<F>>::COSET_SIZE {
                return Err(SynthesisError::InvalidValue(format!("invalid coset size, expected {}, got {}", <I::Combiner as CosetCombiner<F>>::COSET_SIZE, coset_values.len())));
            }
            // let natural_pair_index = (next_domain_idx + (domain_size / 2)) % domain_size;
        
            // for idx in vec![next_domain_idx, natural_pair_index].into_iter() {
            for idx in coset_values.into_iter() {
                let query = iop.query(idx, leaf_values.as_ref());
                queries.push(query);
            }

            roots.push(iop.get_root());

            next_domain_idx = next_domain_idx % (domain_size / 2);
            domain_size >>= 1;
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




// #[test]
// fn test_fib_conversion_into_by_values_fri() {
//     use crate::Fr;
//     use crate::air::Fibonacci;
//     use crate::air::TestTraceSystem;
//     use crate::air::IntoAIR;
//     use crate::arp::*;
//     use crate::ali::ALI;
//     use crate::fft::multicore::Worker;
//     use crate::ali::deep_ali::*;
//     use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;

//     let fib = Fibonacci::<Fr> {
//         final_b: Some(5),
//         at_step: Some(3),
//         _marker: std::marker::PhantomData
//     };

//     let mut test_tracer = TestTraceSystem::<Fr>::new();
//     fib.trace(&mut test_tracer).expect("should work");
//     test_tracer.calculate_witness(1, 1, 3);
//     let mut arp = ARP::<Fr>::new(test_tracer);
//     arp.route_into_single_witness_poly().expect("must work");

//     let mut ali = ALI::from(arp);
//     let alpha = Fr::from_str("123").unwrap();
//     ali.calculate_g(alpha).expect("must work");

//     let mut deep_ali = DeepALI::from(ali);
//     let z = Fr::from_str("62").unwrap();

//     let lde_factor = 8;

//     let worker = Worker::new();

//     println!("F poly size = {}", deep_ali.f_poly.size());
//     println!("G poly size = {}", deep_ali.g_poly.size());

//     let f_lde_values = deep_ali.f_poly.clone().lde(&worker, lde_factor).expect("must work");
//     let g_lde_values = deep_ali.g_poly.clone().lde(&worker, lde_factor).expect("must work");

//     deep_ali.make_deep(f_lde_values, g_lde_values, z).expect("must work");

//     let h1_lde = deep_ali.h_1_poly.take().expect("is something");
//     let h2_lde = deep_ali.h_2_poly.take().expect("is something");

//     let h1_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h1_lde, lde_factor, 1, &worker);
//     let h2_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h2_lde, lde_factor, 1, &worker);

//     // println!("H1 = {:?}", deep_ali.h_1_poly);
//     // println!("H2 = {:?}", deep_ali.h_2_poly);
// }