use ff::PrimeField;
use crate::SynthesisError;
use crate::domains::Domain;
use crate::polynomials::*;

use super::*;
use crate::iop::*;

impl<'a, F: PrimeField, I: IOP<F>> NaiveFriIop<F, I> {
    pub fn verify_prototype(
        proof: & FRIProofPrototype<F, I>,
        leaf_values: & Polynomial<F, Values>, 
        natural_element_index: usize,
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

        let mut expected_value: Option<F> = None;
        let mut domain_size = domain.size as usize;
        let mut next_domain_idx = natural_element_index;

        for (iop_values, iop_challenge) in Some(leaf_values).into_iter().chain(&proof.intermediate_values)
                                        .zip(proof.intermediate_challenges.iter()) {

            let coset_values = <I::Combiner as CosetCombiner<F>>::get_coset_for_natural_index(next_domain_idx, domain_size);

            assert!(coset_values.len() == 2);

            let f_at_omega = I::get_for_natural_index(iop_values.as_ref(), coset_values[0]);

            if let Some(value) = expected_value {
                // check in the next domain
                if *f_at_omega != value {
                    return Ok(false);
                }
            }

            let f_at_minus_omega = I::get_for_natural_index(iop_values.as_ref(), coset_values[1]);
            let divisor = omega_inv.pow([coset_values[0] as u64]);

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
        proof: &FRIProof<F, I>,
        natural_element_index: usize,
        degree: usize, 
        expected_value_from_oracle: F,
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


        let mut expected_value: Option<F> = None;
        let mut domain_size = domain.size as usize;
        let mut next_domain_idx = natural_element_index;

        if proof.queries.len() % degree != 0 {
            return Err(SynthesisError::InvalidValue("invalid number of queries".to_owned()));
        }

        for (round, (root, queries)) in proof.roots.iter()
                                .zip(proof.queries.chunks_exact(degree)) 
                                .enumerate()
        {
            let coset_values = <I::Combiner as CosetCombiner<F>>::get_coset_for_natural_index(next_domain_idx, domain_size);

            if coset_values.len() != <I::Combiner as CosetCombiner<F>>::COSET_SIZE {
                return Err(SynthesisError::InvalidValue(format!("invalid coset size, expected {}, got {}", <I::Combiner as CosetCombiner<F>>::COSET_SIZE, coset_values.len())));
            }

            for q in queries.iter() {
                if !coset_values.contains(&q.natural_index()) {
                    return Ok(false);
                }
            }

            if round == 0 {
                for q in queries.iter() {
                    if q.natural_index() == natural_element_index && q.value() != expected_value_from_oracle {
                        return Ok(false);
                    }
                }
            }

            for (c, q) in coset_values.iter().zip(queries.iter()) {
                let tree_index = <I::Combiner as CosetCombiner<F>>::natural_index_into_tree_index(*c);
                if q.tree_index() != tree_index {
                    return Err(SynthesisError::InvalidValue(format!("invalid tree index for element at natural index {}: expected {}, got {}", c, tree_index, q.tree_index())));
                }
                assert!(q.natural_index() == *c, "coset values and produced queries are expected to be sorted!");
            }

            for q in queries.iter() {
                if !I::verify_query(&q, &root) {
                    return Ok(false);
                }
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
            let divisor = omega_inv.pow([coset_values[0] as u64]);

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
