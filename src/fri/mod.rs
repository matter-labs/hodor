use ff::PrimeField;
use crate::iop::IOP;
use crate::polynomials::*;
use crate::fft::multicore::*;
use crate::SynthesisError;

pub mod fri_on_values;

pub struct FRIIOP<F: PrimeField, I: IOP<F>> {
    _marker_f: std::marker::PhantomData<F>,
    _marker_i: std::marker::PhantomData<I>
}

pub struct FRIProof<F: PrimeField, I: IOP<F>> {
    pub l0_commitment: I,
    pub intermediate_commitments: Vec<I>,
    pub intermediate_values: Vec< Polynomial<F, Values> >,
    pub intermediate_challenges: Vec<F>,
    pub final_coefficients: Vec<F>,
}

impl<F: PrimeField, I: IOP<F>> FRIIOP<F, I> {
    pub fn proof_from_lde(
        lde_values: Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<FRIProof<F, I>, SynthesisError> {
        let l0_commitment: I = I::create(lde_values.as_ref());
        let initial_domain_size = lde_values.size();

        assert!(output_coeffs_at_degree_plus_one.is_power_of_two());
        assert!(lde_factor.is_power_of_two());

        let initial_degree_plus_one = initial_domain_size / lde_factor;
        let num_steps = log2_floor(initial_degree_plus_one / output_coeffs_at_degree_plus_one) as usize;

        let initial_polynomial = lde_values.ifft(&worker);
        let mut initial_polynomial_coeffs = initial_polynomial.into_coeffs();
        initial_polynomial_coeffs.truncate(initial_degree_plus_one);
        
        let mut intermediate_commitments = vec![];
        let mut intermediate_values = vec![];
        let mut intermediate_challenges = vec![];
        let mut next_domain_challenge = l0_commitment.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        let mut next_domain_size = initial_polynomial_coeffs.len() / 2;

        let mut coeffs = initial_polynomial_coeffs;
        
        for _ in 0..num_steps {
            let mut next_coefficients = vec![F::zero(); next_domain_size];
            let coeffs_slice: &[F] = coeffs.as_ref();
            assert!(next_coefficients.len()*2 == coeffs_slice.len());

            worker.scope(next_coefficients.len(), |scope, chunk| {
            for (v, old) in next_coefficients.chunks_mut(chunk)
                            .zip(coeffs_slice.chunks(chunk*2)) {
                scope.spawn(move |_| {
                    for (v, old) in v.iter_mut().zip(old.chunks(2)) {
                        let x = old[0];
                        let mut tmp = old[1];
                        tmp.mul_assign(&next_domain_challenge);
                        tmp.add_assign(&x);

                        *v = tmp;
                    }
                });
            }
        });

        let next_coeffs_as_poly = Polynomial::from_coeffs(next_coefficients.clone())?;
        let next_values_as_poly = next_coeffs_as_poly.lde(&worker, lde_factor)?;
        let intermediate_iop = I::create(next_values_as_poly.as_ref());
        next_domain_challenge = intermediate_iop.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        next_domain_size /= 2;

        intermediate_commitments.push(intermediate_iop);
        intermediate_values.push(next_values_as_poly);

        coeffs = next_coefficients;
    }

    intermediate_challenges.pop().expect("will work");

    assert!(intermediate_challenges.len() == num_steps);
    assert!(intermediate_commitments.len() == num_steps);
    assert!(intermediate_values.len() == num_steps);

    let final_poly_coeffs = coeffs;

    assert!(final_poly_coeffs.len() == output_coeffs_at_degree_plus_one);

    // println!("Final coeffs = {:?}", final_poly_coeffs);

    // let mut degree_plus_one = final_poly_coeffs.len();

    // for v in final_poly_coeffs.iter().rev() {
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

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow+1)) <= num {
        pow += 1;
    }

    pow
}


// #[test]
// fn test_fib_conversion_into_fri() {
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

//     // let h1_coeffs = h1_lde.clone().ifft(&worker);
//     // println!("{:?}", h1_coeffs);

//     let h1_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde(h1_lde.clone(), lde_factor, 1, &worker);
//     let h2_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde(h2_lde.clone(), lde_factor, 1, &worker);

//     // println!("H1 = {:?}", deep_ali.h_1_poly);
//     // println!("H2 = {:?}", deep_ali.h_2_poly);
// }