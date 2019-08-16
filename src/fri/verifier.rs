use ff::PrimeField;
use crate::SynthesisError;
use crate::domains::Domain;
use crate::polynomials::*;

use super::*;
use crate::iop::*;

impl<'a, F: PrimeField, I: IOP<'a, F>> NaiveFriIop<'a, F, I> {
    pub fn verify_prototype(
        proof: &'a FRIProofPrototype<'a, F, I>,
        leaf_values: &'a Polynomial<F, Values>, 
        natural_element_index: usize,
    ) -> Result<bool, SynthesisError>{
        let mut two = F::one();
        two.double();

        let two_inv = two.inverse().ok_or(
            SynthesisError::DivisionByZero(format!("two = {} does not have an inverse", two))
        )?;

        // start from the bottom: we need to get a "pair" and calculate FRI step

        let domain = Domain::<F>::new_for_size((proof.initial_degree_plus_one * proof.lde_factor) as u64)?;

        let domain_element = domain.generator.pow([natural_element_index as u64]);

        let el = domain_element.pow([domain.size]);
        if el != F::one() {
            return Err(SynthesisError::InvalidValue(format!("initial challenge value {} is not in the LDE domain", domain_element)));
        }
        let next_domain_size = domain.size / 2;
        let el = domain_element.pow([next_domain_size]);
        if el == F::one() {
            return Err(SynthesisError::InvalidValue(format!("initial challenge value {} is not in the LDE domain", domain_element)));
        }


        let mut omega = domain.generator;
        let mut omega_inv = omega.inverse().ok_or(
            SynthesisError::DivisionByZero(format!("domain generator {} does not have an inverse", omega))
        )?;

        // let element_x = domain.generator.pow([natural_element_index as u64]);
        // let mut pair_x = element_x;
        // pair_x.negate();

        let mut expected_value: Option<F> = None;
        let mut domain_size = domain.size as usize;
        // let mut next_domain_size = (domain.size / 2) as usize;
        let mut next_domain_idx = natural_element_index;

        for (iop_values, iop_challenge) in Some(leaf_values).into_iter().chain(&proof.intermediate_values)
                                        .zip(proof.intermediate_challenges.iter()) {
            let combiner = I::combine(iop_values.as_ref());

            let natural_pair_index = (next_domain_idx + (domain_size / 2)) % domain_size;

            let f_at_omega = combiner.get(natural_element_index);

            if let Some(value) = expected_value {
                // check in the next domain
                if *f_at_omega != value {
                    return Ok(false);
                }
            }

            let f_at_minus_omega = combiner.get(natural_pair_index as usize);
            let divisor = omega_inv.pow([next_domain_idx as u64]);

            let mut v_even_coeffs = *f_at_omega;
            v_even_coeffs.add_assign(&f_at_minus_omega);

            let mut v_odd_coeffs = *f_at_omega;
            v_odd_coeffs.sub_assign(&f_at_minus_omega);
            v_odd_coeffs.mul_assign(&divisor);

            // those can be treated as (doubled) evaluations of polynomials that
            // are themselves made only from even or odd coefficients of original poly 
            // (with reduction of degree by 2) on a domain of the size twice smaller
            // with an extra factor of "omega" in odd coefficients

            // to do assemble FRI step we just need to add them with a random challenge

            let mut tmp = v_odd_coeffs;
            tmp.mul_assign(&iop_challenge);
            tmp.add_assign(&v_even_coeffs);
            tmp.mul_assign(&two_inv);

            expected_value = Some(tmp);

            next_domain_idx = next_domain_idx % (domain_size / 2);
            domain_size >>= 1;

            omega.square();
            omega_inv.square();
        }


        // finally we need to get expected value from coefficients

        let mut expected_value_from_coefficients = F::zero();
        let mut power = F::one();

        for c in proof.final_coefficients.iter() {
            let mut tmp = power;
            tmp.mul_assign(c);

            expected_value_from_coefficients.add_assign(&tmp);
            power.mul_assign(&omega);
        }
        
        let expected_value = expected_value.expect("is some");

        Ok(expected_value_from_coefficients == expected_value)
    }

    pub fn verify_proof_queries(
        proof: &FRIProof<'a, F, I>,
        natural_element_index: usize,
        degree: usize, 
        expected_value: F,
    ) -> Result<bool, SynthesisError> {
        let mut two = F::one();
        two.double();

        let two_inv = two.inverse().ok_or(
            SynthesisError::DivisionByZero(format!("two = {} does not have an inverse", two))
        )?;

        // start from the bottom: we need to get a "pair" and calculate FRI step

        let domain = Domain::<F>::new_for_size((proof.initial_degree_plus_one * proof.lde_factor) as u64)?;

        let domain_element = domain.generator.pow([natural_element_index as u64]);

        let el = domain_element.pow([domain.size]);
        if el != F::one() {
            return Err(SynthesisError::InvalidValue(format!("initial challenge value {} is not in the LDE domain", domain_element)));
        }
        let next_domain_size = domain.size / 2;
        let el = domain_element.pow([next_domain_size]);
        if el == F::one() {
            return Err(SynthesisError::InvalidValue(format!("initial challenge value {} is not in the LDE domain", domain_element)));
        }

        let mut omega = domain.generator;
        let mut omega_inv = omega.inverse().ok_or(
            SynthesisError::DivisionByZero(format!("domain generator {} does not have an inverse", omega))
        )?;

        let mut expected_value: Option<F> = Some(expected_value);
        let mut domain_size = domain.size as usize;
        let mut next_domain_idx = natural_element_index;

        if proof.queries.len() % degree != 0 {
            return Err(SynthesisError::InvalidValue("invalid number of queries".to_owned()));
        }

        for (root, queries) in proof.roots.iter()
                                .zip(proof.queries.chunks_exact(degree)) 
        {
            if !I::verify_query(&queries[0], &root) {
                return Ok(false);
            }

            if !I::verify_query(&queries[1], &root) {
                return Ok(false);
            }

            if (&queries[0]).index() != next_domain_idx {
                return Ok(false);
            }

            let natural_pair_index = (next_domain_idx + (domain_size / 2)) % domain_size;
            if (&queries[1]).index() != natural_pair_index {
                return Ok(false);
            }

            let iop_challenge = I::encode_root_into_challenge(root);

            let f_at_omega = (&queries[0]).value();
            if let Some(value) = expected_value {
                // check in the next domain
                if f_at_omega != value {
                    return Ok(false);
                }
            }
            let f_at_minus_omega = (&queries[1]).value();
            let divisor = omega_inv.pow([next_domain_idx as u64]);

            let mut v_even_coeffs = f_at_omega;
            v_even_coeffs.add_assign(&f_at_minus_omega);

            let mut v_odd_coeffs = f_at_omega;
            v_odd_coeffs.sub_assign(&f_at_minus_omega);
            v_odd_coeffs.mul_assign(&divisor);

            // those can be treated as (doubled) evaluations of polynomials that
            // are themselves made only from even or odd coefficients of original poly 
            // (with reduction of degree by 2) on a domain of the size twice smaller
            // with an extra factor of "omega" in odd coefficients

            // to do assemble FRI step we just need to add them with a random challenge

            let mut tmp = v_odd_coeffs;
            tmp.mul_assign(&iop_challenge);
            tmp.add_assign(&v_even_coeffs);
            tmp.mul_assign(&two_inv);

            expected_value = Some(tmp);

            next_domain_idx = next_domain_idx % (domain_size / 2);
            domain_size >>= 1;

            omega.square();
            omega_inv.square();
        }


        // finally we need to get expected value from coefficients

        let mut expected_value_from_coefficients = F::zero();
        let mut power = F::one();

        for c in proof.final_coefficients.iter() {
            let mut tmp = power;
            tmp.mul_assign(c);

            expected_value_from_coefficients.add_assign(&tmp);
            power.mul_assign(&omega);
        }
        
        let expected_value = expected_value.expect("is some");

        Ok(expected_value_from_coefficients == expected_value)
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