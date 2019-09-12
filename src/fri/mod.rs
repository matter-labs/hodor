use ff::PrimeField;
use crate::iop::*;
use crate::polynomials::*;
use crate::fft::multicore::*;
use crate::SynthesisError;
use crate::utils::log2_floor;

pub mod fri_on_values;
pub mod verifier;
pub mod query_producer;

pub trait FriProofPrototype<F: PrimeField, I: IOP<F>> {
    fn get_roots(&self) -> Vec< < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput>;
    fn get_final_root(&self) -> < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput;
    fn get_final_coefficients(&self) -> Vec<F>;
}

pub trait FriProof<F: PrimeField, I: IOP<F>> {
    fn get_final_coefficients(&self) -> &[F];
}

pub trait FriIop<F: PrimeField> {
    const DEGREE: usize;

    type IopType: IOP<F>;
    type ProofPrototype: FriProofPrototype<F, Self::IopType>;
    type Proof: FriProof<F, Self::IopType>;

    fn proof_from_lde(
        lde_values: &Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<Self::ProofPrototype, SynthesisError>;

    fn prototype_into_proof(
        prototype: Self::ProofPrototype,
        iop_values: &Polynomial<F, Values>,
        natural_first_element_index: usize,
    ) -> Result<Self::Proof, SynthesisError>;

    fn verify_proof(
        proof: &Self::Proof,
        natural_element_index: usize,
        expected_value: F
    ) -> Result<bool, SynthesisError>;
}

pub struct NaiveFriIop<F: PrimeField, I: IOP<F>> {
    _marker_f: std::marker::PhantomData<F>,
    _marker_i: std::marker::PhantomData<I>
}

impl<'a, F: PrimeField, I: IOP<F>> FriIop<F> for NaiveFriIop<F, I> {
    const DEGREE: usize = 2;

    type IopType = I;
    type ProofPrototype = FRIProofPrototype<F, I>;
    type Proof = FRIProof<F, I>;

    fn proof_from_lde(
        lde_values: &Polynomial<F, Values>, 
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
        iop_values: &Polynomial<F, Values>,
        natural_first_element_index: usize,
    ) -> Result<Self::Proof, SynthesisError> {
        prototype.produce_proof(iop_values, natural_first_element_index)
    }

    fn verify_proof(
        proof: &Self::Proof,
        natural_element_index: usize,
        expected_value: F
    ) -> Result<bool, SynthesisError> {
        Self::verify_proof_queries(proof, natural_element_index, Self::DEGREE, expected_value)
    }
}

#[derive(PartialEq, Eq, Clone)]
pub struct FRIProofPrototype<F: PrimeField, I: IOP<F>> {
    pub l0_commitment: I,
    pub intermediate_commitments: Vec<I>,
    pub intermediate_values: Vec< Polynomial<F, Values> >,
    pub challenges: Vec<F>,
    pub final_root: < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput,
    pub final_coefficients: Vec<F>,
    pub initial_degree_plus_one: usize,
    pub output_coeffs_at_degree_plus_one: usize,
    pub lde_factor: usize,
}

impl<F: PrimeField, I: IOP<F>> FriProofPrototype<F, I> for FRIProofPrototype<F, I> {
    fn get_roots(&self) -> Vec< < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput> {
        let mut roots = vec![];
        roots.push(self.l0_commitment.get_root().clone());
        for c in self.intermediate_commitments.iter() {
            roots.push(c.get_root().clone());
        }

        roots
    }

    fn get_final_root(&self) -> < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput {
        self.final_root.clone()
    }

    fn get_final_coefficients(&self) -> Vec<F> {
        self.final_coefficients.clone()
    }
}

#[derive(PartialEq, Eq, Clone)]
pub struct FRIProof<F: PrimeField, I: IOP<F>> {
    pub queries: Vec< <I as IOP<F> >::Query >,
    pub roots: Vec< < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput>,
    pub final_coefficients: Vec<F>,
    pub initial_degree_plus_one: usize,
    pub output_coeffs_at_degree_plus_one: usize,
    pub lde_factor: usize,
}

impl<F: PrimeField, I: IOP<F>> FriProof<F, I> for FRIProof<F, I> {
    fn get_final_coefficients(&self) -> &[F] {
        &self.final_coefficients
    }
}

impl<'a, F: PrimeField, I: IOP<F>> NaiveFriIop<F, I> {
    pub fn proof_from_lde_through_coefficients(
        lde_values: Polynomial<F, Values>, 
        lde_factor: usize,
        output_coeffs_at_degree_plus_one: usize,
        worker: &Worker
    ) -> Result<FRIProofPrototype<F, I>, SynthesisError> {
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
        let mut challenges = vec![];
        let mut next_domain_challenge = l0_commitment.get_challenge_scalar_from_root();
        challenges.push(next_domain_challenge);
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
                        // a_0 + beta * a_1

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
        challenges.push(next_domain_challenge);

        next_domain_size /= 2;

        intermediate_commitments.push(intermediate_iop);
        intermediate_values.push(next_values_as_poly);

        coeffs = next_coefficients;
    }

    challenges.pop().expect("will work");

    let final_root = roots.pop().expect("will work");

    assert!(challenges.len() == num_steps);
    assert!(intermediate_commitments.len() == num_steps);
    assert!(intermediate_values.len() == num_steps);

    let final_poly_coeffs = coeffs;

    assert!(final_poly_coeffs.len() == output_coeffs_at_degree_plus_one);

    Ok(FRIProofPrototype {
        l0_commitment,
        intermediate_commitments,
        intermediate_values,
        challenges,
        final_root,
        final_coefficients: final_poly_coeffs,
        initial_degree_plus_one,
        output_coeffs_at_degree_plus_one,
        lde_factor,
    })

    }
}

#[test]
fn test_one_fri_step() {
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

    let mut lde_coeffs = vec![];
    let mut f = Fr::one();
    for _ in 0..4 {
        lde_coeffs.push(f);
        f.double();
    }

    let lde_factor = 4;
    let worker = Worker::new();

    let lde_coeffs = Polynomial::<Fr, Coefficients>::from_coeffs(lde_coeffs).expect("must work");

    let lde_values = lde_coeffs.clone().lde(&worker, lde_factor).unwrap();

    let output_at_degree_plus_one = 2;

    let fri_proof_by_values = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&lde_values, lde_factor, output_at_degree_plus_one, &worker).expect("must work");
    let fri_proof_by_coeffs = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_through_coefficients(lde_values.clone(), lde_factor, output_at_degree_plus_one, &worker).expect("must work");
        
    let coset_index = 3;
    let coset_pair_index = coset_index + lde_factor * 2;

    let generator_inv = lde_values.omegainv;
    let divisor = generator_inv.pow([coset_index as u64]);

    let two_inv = Fr::from_str("2").unwrap().inverse().unwrap();

    let challenge = fri_proof_by_values.challenges[0];

    let value_at_omega = lde_values.as_ref()[coset_index];
    let value_at_minus_omega = lde_values.as_ref()[coset_pair_index];

    // interpolation

    let mut t0 = value_at_omega;
    t0.add_assign(&value_at_minus_omega);

    let mut t1 = value_at_omega;
    t1.sub_assign(&value_at_minus_omega);
    t1.mul_assign(&divisor);

    t1.mul_assign(&challenge);

    t0.add_assign(&t1);
    t0.mul_assign(&two_inv);

    let mut new_coeffs = vec![];
    for c in lde_coeffs.as_ref().chunks(2) {
        let mut tmp = c[1];
        tmp.mul_assign(&challenge);
        tmp.add_assign(&c[0]);
        new_coeffs.push(tmp);
    }

    assert!(fri_proof_by_values.final_coefficients == new_coeffs);

    let next_lde = Polynomial::from_coeffs(new_coeffs).unwrap().lde(&worker, lde_factor).unwrap();

    let next_domain_index = coset_index;

    let value = next_lde.as_ref()[next_domain_index];

    assert!(value == t0);

    let lde = lde_values;

    for i in (1..lde.size()).step_by(2) {


        assert!(fri_proof_by_values.final_coefficients == fri_proof_by_coeffs.final_coefficients);
        assert!(fri_proof_by_values.final_root == fri_proof_by_coeffs.final_root);
        assert!(fri_proof_by_values.intermediate_values == fri_proof_by_coeffs.intermediate_values);
        assert!(fri_proof_by_values.challenges == fri_proof_by_coeffs.challenges);
        assert!(fri_proof_by_values.l0_commitment == fri_proof_by_coeffs.l0_commitment);
        assert!(fri_proof_by_values.intermediate_commitments == fri_proof_by_coeffs.intermediate_commitments);
            
        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &fri_proof_by_coeffs,
            &lde,
            i
        );

        if valid.is_err() {
            println!("failed for i = {}", i);
            println!("error = {}", valid.err().unwrap())
        } else {
            let valid = valid.unwrap();
            if !valid {
                println!("failed for i = {}", i);
            }
        }
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

    let natural_x_index = 63;
    // let natural_x_index = 33;
    // let natural_x_index = 1;

    {
        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &h1_fri_proof,
            &h1_lde,
            natural_x_index
        ).expect("must work");

        assert!(valid);
    }

    {
        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &h2_fri_proof,
            &h2_lde,
            natural_x_index
        ).expect("must work");

        assert!(valid);
    }

    println!("Making queries and true proofs");

    let proof_h1 = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr> >::prototype_into_proof(
        h1_fri_proof,
        &h1_lde,
        natural_x_index
    ).expect("must work");

    let expected_value = h1_lde.as_ref()[natural_x_index];

    let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr> >::verify_proof(
        &proof_h1,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(valid);

    let proof_h2 = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FriIop<Fr>>::prototype_into_proof(
        h2_fri_proof,
        &h2_lde,
        natural_x_index
    ).expect("must work");

    let expected_value = h2_lde.as_ref()[natural_x_index];

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FriIop<Fr>>::verify_proof(
        &proof_h2,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(valid);

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FriIop<Fr>>::verify_proof(
        &proof_h1,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(!valid);

}

#[test]
fn test_fri_on_values_vs_on_coefficients() {
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

    for i in (17..h1_lde.size()).step_by(2) {

        // let h1_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h1_lde, lde_factor, output_at_degree, &worker).expect("must work");
        let h1_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_through_coefficients(h1_lde.clone(), lde_factor, output_at_degree, &worker).expect("must work");
        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &h1_fri_proof,
            &h1_lde,
            i
        );

        if valid.is_err() {
            println!("failed for i = {}", i);
            println!("error = {}", valid.err().unwrap())
        } else {
            let valid = valid.unwrap();
            if !valid {
                println!("failed for i = {}", i);
            }
        }
    }

    let h1_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h1_lde, lde_factor, output_at_degree, &worker).expect("must work");
    let h2_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h2_lde, lde_factor, output_at_degree, &worker).expect("must work");

    let h1_fri_proof_from_coeffs = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_through_coefficients(h1_lde.clone(), lde_factor, output_at_degree, &worker).expect("must work");
    let h2_fri_proof_from_coeffs = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_through_coefficients(h2_lde.clone(), lde_factor, output_at_degree, &worker).expect("must work");

    assert!(h1_fri_proof_from_coeffs.final_coefficients == h1_fri_proof.final_coefficients);
    assert!(h1_fri_proof_from_coeffs.final_root == h1_fri_proof.final_root);
    assert!(h1_fri_proof_from_coeffs.intermediate_values == h1_fri_proof.intermediate_values);
    assert!(h1_fri_proof_from_coeffs.challenges == h1_fri_proof.challenges);
    assert!(h1_fri_proof_from_coeffs.l0_commitment == h1_fri_proof.l0_commitment);
    assert!(h1_fri_proof_from_coeffs.intermediate_commitments == h1_fri_proof.intermediate_commitments);

    let natural_x_index = 63;
    // let natural_x_index = 33;
    // let natural_x_index = 1;

    {
        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &h1_fri_proof_from_coeffs,
            &h1_lde,
            natural_x_index
        ).expect("must work");

        assert!(valid);

        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &h1_fri_proof,
            &h1_lde,
            natural_x_index
        ).expect("must work");

        assert!(valid);
    }

    {
        let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::verify_prototype(
            &h2_fri_proof,
            &h2_lde,
            natural_x_index
        ).expect("must work");

        assert!(valid);
    }

    println!("Making queries and true proofs");

    let proof_h1 = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr> >::prototype_into_proof(
        h1_fri_proof,
        &h1_lde,
        natural_x_index
    ).expect("must work");

    let expected_value = h1_lde.as_ref()[natural_x_index];

    let valid = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr> >::verify_proof(
        &proof_h1,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(valid);

    let proof_h2 = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FriIop<Fr>>::prototype_into_proof(
        h2_fri_proof,
        &h2_lde,
        natural_x_index
    ).expect("must work");

    let expected_value = h2_lde.as_ref()[natural_x_index];

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FriIop<Fr>>::verify_proof(
        &proof_h2,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(valid);

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FriIop<Fr>>::verify_proof(
        &proof_h1,
        natural_x_index,
        expected_value
    ).expect("must work");

    assert!(!valid);

}