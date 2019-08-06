use ff::PrimeField;
use crate::iop::IOP;
use crate::polynomials::*;
use crate::fft::multicore::*;
use crate::SynthesisError;
use crate::utils::*;

use super::{FRIIOP, FRIProof};

impl<F: PrimeField, I: IOP<F>> FRIIOP<F, I> {
    pub fn proof_from_lde_by_values(
        lde_values: &Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<FRIProof<F, I>, SynthesisError> {
        let l0_commitment: I = I::create(lde_values.as_ref());
        let initial_domain_size = lde_values.size();
        let precomputation_size = initial_domain_size/2;
        let mut two = F::one();
        two.double();
        let two_inv = two.inverse().expect("should exist");

        let mut omegas_inv = vec![F::zero(); precomputation_size];
        let omega_inv = lde_values.omegainv;

        // for interpolations we will need factors 2*w^k in denominator,
        // so we just make powers

        worker.scope(omegas_inv.len(), |scope, chunk| {
            for (i, v) in omegas_inv.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = omega_inv.pow(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        *v = u;
                        u.mul_assign(&omega_inv);
                    }
                });
            }
        });

        assert!(output_coeffs_at_degree_plus_one.is_power_of_two());
        assert!(lde_factor.is_power_of_two());

        let initial_degree_plus_one = initial_domain_size / lde_factor;
        let num_steps = log2_floor(initial_degree_plus_one / output_coeffs_at_degree_plus_one) as usize;

        let mut intermediate_commitments = vec![];
        let mut intermediate_values = vec![];
        let mut intermediate_challenges = vec![];
        let mut next_domain_challenge = l0_commitment.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        let mut next_domain_size = initial_domain_size / 2;
        let mut stride = 1;

        let mut values_slice = lde_values.as_ref();

        let omegas_inv_ref: &[F] = omegas_inv.as_ref();
        
        for _ in 0..num_steps {
            let mut next_values = vec![F::zero(); next_domain_size];

            assert!(values_slice.len() == next_values.len() * 2);

            worker.scope(next_values.len(), |scope, chunk| {
            for (i, v) in next_values.chunks_mut(chunk)
                            .enumerate() {
                scope.spawn(move |_| {
                    let initial_k = i*chunk;
                    for (j, v) in v.iter_mut().enumerate() {
                        let idx = initial_k + j;
                        let omega_idx = idx * stride;
                        let f_at_omega = values_slice[idx];
                        let f_at_minus_omega = values_slice[idx + next_domain_size];

                        let mut v_even_coeffs = f_at_omega;
                        v_even_coeffs.add_assign(&f_at_minus_omega);

                        let mut v_odd_coeffs = f_at_omega;
                        v_odd_coeffs.sub_assign(&f_at_minus_omega);
                        v_odd_coeffs.mul_assign(&omegas_inv_ref[omega_idx]);

                        // those can be treated as (doubled) evaluations of polynomials that
                        // are themselves made only from even or odd coefficients of original poly 
                        // (with reduction of degree by 2) on a domain of the size twice smaller
                        // with an extra factor of "omega" in odd coefficients

                        // to do assemble FRI step we just need to add them with a random challenge

                        let mut tmp = v_odd_coeffs;
                        tmp.mul_assign(&next_domain_challenge);
                        tmp.add_assign(&v_even_coeffs);
                        tmp.mul_assign(&two_inv);

                        *v = tmp;
                    }
                });
            }
        });

        let intermediate_iop = I::create(next_values.as_ref());
        next_domain_challenge = intermediate_iop.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        next_domain_size /= 2;
        stride *= 2;

        intermediate_commitments.push(intermediate_iop);
        let next_values_as_poly = Polynomial::from_values(next_values)?;
        intermediate_values.push(next_values_as_poly);

        values_slice = intermediate_values.last().expect("is something").as_ref();
    }

    intermediate_challenges.pop().expect("will work");

    assert!(intermediate_challenges.len() == num_steps);
    assert!(intermediate_commitments.len() == num_steps);
    assert!(intermediate_values.len() == num_steps);

    let final_poly_values = Polynomial::from_values(values_slice.to_vec())?;
    let final_poly_coeffs = final_poly_values.ifft(&worker);

    let mut final_poly_coeffs = final_poly_coeffs.into_coeffs();
    final_poly_coeffs.truncate(output_coeffs_at_degree_plus_one);

    // println!("Final coeffs = {:?}", final_poly_coeffs.as_ref());

    // let mut degree_plus_one = final_poly_coeffs.size();

    // for v in final_poly_coeffs.as_ref().iter().rev() {
    //     if v.is_zero() {
    //         degree_plus_one -= 1;
    //     } else {
    //         break;
    //     }
    // }

    // println!("Degree = {}", degree_plus_one);

    Ok(FRIProof {
        l0_commitment,
        intermediate_commitments,
        intermediate_values,
        intermediate_challenges,
        final_coefficients: final_poly_coeffs
    })

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