use ff::PrimeField;
use crate::iop::IOP;
use crate::polynomials::*;
use crate::fft::multicore::*;
use crate::SynthesisError;

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

        let mut omegas = vec![F::zero(); precomputation_size];
        let omega = lde_values.omega;

        let mut two_omegas_inv = vec![F::zero(); precomputation_size];
        let omega_inv = lde_values.omegainv;

        // for interpolations we will need factors 2*w^k in denominator,
        // so we just make powers

        worker.scope(two_omegas_inv.len(), |scope, chunk| {
            for (i, v) in two_omegas_inv.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = omega_inv.pow(&[(i * chunk) as u64]);
                    u.mul_assign(&two_inv);
                    for v in v.iter_mut() {
                        *v = u;
                        u.mul_assign(&omega_inv);
                    }
                });
            }
        });

        worker.scope(omegas.len(), |scope, chunk| {
            for (i, v) in omegas.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = omega.pow(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        *v = u;
                        u.mul_assign(&omega);
                    }
                });
            }
        });

        assert!(output_coeffs_at_degree_plus_one.is_power_of_two());
        assert!(lde_factor.is_power_of_two());

        let initial_degree_plus_one = initial_domain_size / lde_factor;
        println!("Original poly degree = {}, will output at degree = {}", initial_degree_plus_one, output_coeffs_at_degree_plus_one);
        let num_steps = log2_floor(initial_degree_plus_one / output_coeffs_at_degree_plus_one) as usize;
        println!("Will make {} FRI steps", num_steps);

        let mut intermediate_commitments = vec![];
        let mut intermediate_values = vec![];
        let mut intermediate_challenges = vec![];
        let mut next_domain_challenge = l0_commitment.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        let mut next_domain_size = initial_domain_size / 2;
        let mut stride = 1;

        let mut values_slice = lde_values.as_ref();
        let omegas_ref: &[F] = omegas.as_ref();
        let two_omegas_inv_ref: &[F] = two_omegas_inv.as_ref();
        
        for _ in 0..num_steps {
            let mut next_values = vec![F::zero(); next_domain_size];

            worker.scope(next_values.len(), |scope, chunk| {
            for (i, v) in next_values.chunks_mut(chunk)
                            .enumerate() {
                scope.spawn(move |_| {
                    let initial_k = i*chunk;
                    // use naive linear interpolation like
                    // (f(omega^k) - f(omega^(k + N/2)))/(omega^k - omega^(k + N/2))*(Y - omega^(k + N/2)) + f(omega^(k + N/2)) =
                    // (f(omega^k) - f(omega^(k + N/2)))/(2*omega^k)*(Y + omega^(k)) + f(omega^(k + N/2))
                    // and immediately evaluate
                    for (j, v) in v.iter_mut().enumerate() {
                        // f(omega^k)
                        let mut alpha = values_slice[initial_k + j];
                        // f(omega^k) - f(omega^(k + N/2))
                        alpha.sub_assign(&values_slice[initial_k + j + next_domain_size]);
                        // (f(omega^k) - f(omega^(k + N/2)))/(2*omega^k)
                        alpha.mul_assign(&two_omegas_inv_ref[initial_k + j*stride]);

                        // Y
                        let mut y = next_domain_challenge;
                        // (Y + omega^(k))
                        y.add_assign(&omegas_ref[initial_k + j*stride]);

                        // alpha * (Y + omega^(k))
                        alpha.mul_assign(&y);
                        // + f(omega^(k + N/2))
                        alpha.add_assign(&values_slice[initial_k + j + next_domain_size]);
                        *v = alpha;

                        // // check the interpolation
                        // {
                        //     let c = values_slice[initial_k + j + next_domain_size];
                        //     let mut x = omegas_ref[initial_k + j*stride];
                        //     x.double();
                        //     let mut t = alpha;
                        //     t.mul_assign(&x);
                        //     t.add_assign(&c);
                        //     debug_assert!(t == values_slice[initial_k + j]);


                        // }
                    }
                });
            }
        });

        let intermediate_iop = I::create(next_values.as_ref());
        next_domain_challenge = intermediate_iop.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        next_domain_size /= 2;
        stride += 1;

        intermediate_commitments.push(intermediate_iop);
        let next_values_as_poly = Polynomial::from_values(next_values)?;
        let next_values_as_coeffs = next_values_as_poly.clone().ifft(&worker);
        println!("Next poly = {:?}", next_values_as_coeffs.as_ref());
        intermediate_values.push(next_values_as_poly);

        values_slice = intermediate_values.last().expect("is something").as_ref();
        }

    intermediate_challenges.pop().expect("will work");

    assert!(intermediate_challenges.len() == num_steps);
    assert!(intermediate_commitments.len() == num_steps);
    assert!(intermediate_values.len() == num_steps);

    let final_poly_values = Polynomial::from_values(values_slice.to_vec())?;
    let final_poly_coeffs = final_poly_values.ifft(&worker);

    println!("Final coeffs = {:?}", final_poly_coeffs.as_ref());

    let mut degree_plus_one = final_poly_coeffs.size();

    for v in final_poly_coeffs.as_ref().iter().rev() {
        if v.is_zero() {
            degree_plus_one -= 1;
        } else {
            break;
        }
    }

    println!("Degree = {}", degree_plus_one);

    Ok(FRIProof {
        l0_commitment,
        intermediate_commitments,
        intermediate_values,
        intermediate_challenges,
        final_coefficients: final_poly_coeffs.into_coeffs()
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


#[test]
fn test_fib_conversion_into_fri() {
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::arp::*;
    use crate::ali::ALI;
    use crate::fft::multicore::Worker;
    use crate::ali::deep_ali::*;
    use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let mut arp = ARP::<Fr>::new(test_tracer);
    arp.route_into_single_witness_poly().expect("must work");

    let mut ali = ALI::from(arp);
    let alpha = Fr::from_str("123").unwrap();
    ali.calculate_g(alpha).expect("must work");

    let mut deep_ali = DeepALI::from(ali);
    let z = Fr::from_str("62").unwrap();

    let lde_factor = 8;

    let worker = Worker::new();

    println!("F poly size = {}", deep_ali.f_poly.size());
    println!("G poly size = {}", deep_ali.g_poly.size());

    let f_lde_values = deep_ali.f_poly.clone().lde(&worker, lde_factor).expect("must work");
    let g_lde_values = deep_ali.g_poly.clone().lde(&worker, lde_factor).expect("must work");

    deep_ali.make_deep(f_lde_values, g_lde_values, z).expect("must work");

    let h1_lde = deep_ali.h_1_poly.take().expect("is something");
    let h2_lde = deep_ali.h_2_poly.take().expect("is something");

    let h1_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde(&h1_lde, lde_factor, 1, &worker);
    let h2_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde(&h2_lde, lde_factor, 1, &worker);

    // println!("H1 = {:?}", deep_ali.h_1_poly);
    // println!("H2 = {:?}", deep_ali.h_2_poly);
}