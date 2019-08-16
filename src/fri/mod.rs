use ff::PrimeField;
use crate::iop::*;
use crate::polynomials::*;
use crate::fft::multicore::*;
use crate::SynthesisError;
use crate::utils::log2_floor;

pub mod fri_on_values;
pub mod verifier;
pub mod query_producer;

pub trait FRIIOP<'a, F: PrimeField> {
    const DEGREE: usize;

    type IopType: IOP<'a, F>;
    type ProofPrototype;
    type Proof;

    fn proof_from_lde(
        lde_values: &'a Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<Self::ProofPrototype, SynthesisError>;

    fn prototype_into_proof(
        prototype: Self::ProofPrototype,
        iop_values: &'a Polynomial<F, Values>,
        natural_first_element_index: usize,
    ) -> Result<Self::Proof, SynthesisError>;

    fn verify_proof(
        proof: &'a Self::Proof,
        natural_element_index: usize,
        expected_value: F
    ) -> Result<bool, SynthesisError>;
}

pub struct NaiveFriIop<'a, F: PrimeField, I: IOP<'a, F>> {
    _marker_f: std::marker::PhantomData<F>,
    _marker_i: std::marker::PhantomData<&'a I>
}

impl<'a, F: PrimeField, I: IOP<'a, F>> FRIIOP<'a, F> for NaiveFriIop<'a, F, I> {
    const DEGREE: usize = 2;

    type IopType = I;
    type ProofPrototype = FRIProofPrototype<'a, F, I>;
    type Proof = FRIProof<'a, F, I>;

    fn proof_from_lde(
        lde_values: &'a Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<Self::ProofPrototype, SynthesisError> {
        NaiveFriIop::proof_from_lde_by_values(
            lde_values, 
            lde_factor,
            output_coeffs_at_degree_plus_one,
            worker
        )
    }

    fn prototype_into_proof(
        prototype: Self::ProofPrototype,
        iop_values: &'a Polynomial<F, Values>,
        natural_first_element_index: usize,
    ) -> Result<Self::Proof, SynthesisError> {
        prototype.produce_proof(iop_values, natural_first_element_index)
    }

    fn verify_proof(
        proof: &'a Self::Proof,
        natural_element_index: usize,
        expected_value: F
    ) -> Result<bool, SynthesisError> {
        Self::verify_proof_queries(proof, natural_element_index, Self::DEGREE, expected_value)
    }
}

pub struct FRIProofPrototype<'a, F: PrimeField, I: IOP<'a, F>> {
    pub l0_commitment: I,
    pub intermediate_commitments: Vec<I>,
    pub intermediate_values: Vec< Polynomial<F, Values> >,
    pub intermediate_challenges: Vec<F>,
    pub final_root: < <I::Tree as IopTree<'a, F> >::Hasher as IopTreeHasher<F>>::HashOutput,
    pub final_coefficients: Vec<F>,
    pub initial_degree_plus_one: usize,
    pub output_coeffs_at_degree_plus_one: usize,
    pub lde_factor: usize,
    _marker: std::marker::PhantomData<&'a I>
}

pub struct FRIProof<'a, F: PrimeField, I: IOP<'a, F>> {
    pub queries: Vec< <I as IOP<'a, F> >::Query >,
    pub roots: Vec< < <I::Tree as IopTree<'a, F> >::Hasher as IopTreeHasher<F>>::HashOutput>,
    pub final_coefficients: Vec<F>,
    pub initial_degree_plus_one: usize,
    pub output_coeffs_at_degree_plus_one: usize,
    pub lde_factor: usize,
    _marker: std::marker::PhantomData<&'a I>
}

impl<'a, F: PrimeField, I: IOP<'a, F>> NaiveFriIop<'a, F, I> {
    pub fn proof_from_lde_through_coefficients(
        lde_values: Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<FRIProofPrototype<'a, F, I>, SynthesisError> {
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
        let mut roots = vec![];
        
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
        let root = intermediate_iop.get_root();
        roots.push(root);
        next_domain_challenge = intermediate_iop.get_challenge_scalar_from_root();
        intermediate_challenges.push(next_domain_challenge);
        next_domain_size /= 2;

        intermediate_commitments.push(intermediate_iop);
        intermediate_values.push(next_values_as_poly);

        coeffs = next_coefficients;
    }

    intermediate_challenges.pop().expect("will work");

    let final_root = roots.pop().expect("will work");

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

    Ok(FRIProofPrototype {
        l0_commitment,
        intermediate_commitments,
        intermediate_values,
        intermediate_challenges,
        final_root,
        final_coefficients: final_poly_coeffs,
        initial_degree_plus_one,
        output_coeffs_at_degree_plus_one,
        lde_factor,
        _marker: std::marker::PhantomData
    })

    }
}

#[test]
fn test_fib_fri_iop_verifier() {
    use ff::Field;
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::fft::multicore::Worker;
    use crate::transcript::Transcript;
    use crate::arp::*;
    use crate::ali::*;
    use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;
    use crate::iop::blake2s_trivial_iop::Blake2sIopTree;
    use crate::fri::*;
    use crate::transcript::*;
    use crate::ali::per_register::*;
    use crate::iop::*;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let lde_factor = 16;
    let mut transcript = Blake2sTranscript::new();

    let worker = Worker::new();

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();
    let witness = witness.expect("some witness");
    // println!("Witness = {:?}", witness);

    let is_satisfied = ARPInstance::<Fr, PerRegisterARP>::is_satisfied(&props, &witness, &worker);
    assert!(is_satisfied.is_ok());

    let arp = ARPInstance::<Fr, PerRegisterARP>::from_instance(props.clone(), &worker).expect("must work");

    let witness_polys = arp.calculate_witness_polys(witness, &worker).expect("must work");

    let f_ldes: Vec<_> = witness_polys.iter().map(|w| {
        w.clone().lde(&worker, lde_factor).expect("must work")
    }).collect();

    let f_oracles: Vec<_> = f_ldes.iter().map(|l|
        Blake2sIopTree::create(l.as_ref())
    ).collect(); 

    for o in f_oracles.iter() {
        transcript.commit_bytes(&o.get_root()[..]);
    }

    let ali = ALIInstance::from_arp(arp, &worker).expect("is some");

    let g_poly_interpolant = ali.calculate_g(&mut transcript, witness_polys.clone(), &worker).expect("is some");

    let g_lde = g_poly_interpolant.clone().lde(&worker, lde_factor).expect("is something");

    let g_oracle = Blake2sIopTree::create(g_lde.as_ref());
    transcript.commit_bytes(&g_oracle.get_root());

    let (h1_lde, h2_lde, _f_at_z_m, _g_at_z) = ali.calculate_deep(
        &witness_polys,
        &f_ldes,
        &g_poly_interpolant,
        &g_lde,
        &mut transcript,
        &worker
    ).expect("must work");

    let output_at_degree = 1;

    let h1_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h1_lde, lde_factor, output_at_degree, &worker).expect("must work");
    let h2_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h2_lde, lde_factor, output_at_degree, &worker).expect("must work");

    let natural_x_index = 1;

    // {
    //     let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
    //         &h1_fri_proof,
    //         &h1_lde,
    //         natural_x_index
    //     ).expect("must work");

    //     assert!(valid);
    // }

    // {
    //     let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
    //         &h2_fri_proof,
    //         &h2_lde,
    //         natural_x_index
    //     ).expect("must work");

    //     assert!(valid);
    // }

    let proof_h1 = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::prototype_into_proof(
        h1_fri_proof,
        &h1_lde,
        natural_x_index
    ).expect("must work");

    let expected_value = h1_lde.as_ref()[natural_x_index];

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::verify_proof(
        &proof_h1,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(valid);

    let proof_h2 = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::prototype_into_proof(
        h2_fri_proof,
        &h2_lde,
        natural_x_index
    ).expect("must work");

    let expected_value = h2_lde.as_ref()[natural_x_index];

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::verify_proof(
        &proof_h2,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(valid);

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::verify_proof(
        &proof_h1,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(!valid);

}