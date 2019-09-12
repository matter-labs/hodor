
use ff::PrimeField;

use crate::air::IntoAIR;
use crate::fft::multicore::Worker;
use crate::transcript::Transcript;
use crate::arp::*;
use crate::ali::*;
use crate::fri::*;
use crate::transcript::*;
use crate::ali::per_register::*;
use crate::iop::*;
use crate::fri::*;
use crate::SynthesisError;
use crate::verifier::*;

pub struct Prover<F: PrimeField, T: Transcript<F>, I: IOP<F>, P: FriProofPrototype<F, I>, PR: FriProof<F, I>, FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR>, A: ARPType> {
    arp: ARPInstance::<F, A>,
    ali: ALIInstance::<F, A>,
    lde_factor: usize,
    fri_final_degree_plus_one: usize,
    worker: Worker,

    _marker_t: std::marker::PhantomData<T>,
    _marker_i: std::marker::PhantomData<I>,
    _marker_p: std::marker::PhantomData<P>,
    _marker_pr: std::marker::PhantomData<PR>,
    _marker_fri: std::marker::PhantomData<FRI>,

}

impl<F: PrimeField, T: Transcript<F>, I: IOP<F>, P: FriProofPrototype<F, I>, PR: FriProof<F, I>, FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR> > Prover<F, T, I, P, PR, FRI, PerRegisterARP> {

    pub fn new(instance: InstanceProperties<F>, lde_factor: usize, fri_final_degree_plus_one: usize) -> Result<Self, SynthesisError> {
        let worker = Worker::new();
        let arp = ARPInstance::<F, PerRegisterARP>::from_instance(instance, &worker)?;
        let ali = ALIInstance::from_arp(arp.clone(), &worker)?;

        Ok(Self {
            arp,
            ali,
            lde_factor,
            fri_final_degree_plus_one,
            worker,

            _marker_t: std::marker::PhantomData,
            _marker_i: std::marker::PhantomData,
            _marker_p: std::marker::PhantomData,
            _marker_pr: std::marker::PhantomData,
            _marker_fri: std::marker::PhantomData,
        })
    }

    pub fn prove(&self, witness: Vec<Vec<F>>) -> Result<InstanceProof<F, T, I, P, PR, FRI, PerRegisterARP>, SynthesisError> {
        let mut transcript = T::new();

        let witness_polys = self.arp.calculate_witness_polys(witness, &self.worker)?;

        let mut f_ldes = Vec::with_capacity(witness_polys.len());

        for w in witness_polys.iter() {
            let f_lde = w.clone().lde(&self.worker, self.lde_factor)?;
            f_ldes.push(f_lde);
        }

        let f_oracles: Vec<_> = f_ldes.iter().map(|l|
            I::create(l.as_ref())
        ).collect(); 

        let mut f_iop_roots = vec![];
        for o in f_oracles.iter() {
            let root = o.get_root();
            transcript.commit_bytes(root.as_ref());
            f_iop_roots.push(root);
        }

        let g_poly_interpolant = self.ali.calculate_g(&mut transcript, witness_polys.clone(), &self.worker)?;

        let g_lde = g_poly_interpolant.clone().lde(&self.worker, self.lde_factor)?;

        let g_oracle = I::create(g_lde.as_ref());
        let g_iop_root = g_oracle.get_root();
        transcript.commit_bytes(g_iop_root.as_ref());

        println!("Calculating DEEP part");

        let (h1_lde, h2_lde, f_at_z_m, _g_at_z) = self.ali.calculate_deep(
            &witness_polys,
            &f_ldes,
            &g_poly_interpolant,
            &g_lde,
            &mut transcript,
            &self.worker
        )?;

        let fri_final_poly_degree = self.fri_final_degree_plus_one;

        println!("Calculating FRI part");

        let h1_fri_proof_proto = FRI::proof_from_lde(&h1_lde, self.lde_factor, fri_final_poly_degree, &self.worker)?;
        let h2_fri_proof_proto = FRI::proof_from_lde(&h2_lde, self.lde_factor, fri_final_poly_degree, &self.worker)?;

        let h1_iop_roots = h1_fri_proof_proto.get_roots();
        let h2_iop_roots = h2_fri_proof_proto.get_roots();

        // TODO: we can also potentially commit intermediate roots

        transcript.commit_bytes(&h1_fri_proof_proto.get_final_root().as_ref());
        for el in h1_fri_proof_proto.get_final_coefficients().into_iter() {
            transcript.commit_field_element(&el);
        }
        transcript.commit_bytes(&h2_fri_proof_proto.get_final_root().as_ref());
        for el in h2_fri_proof_proto.get_final_coefficients().into_iter() {
            transcript.commit_field_element(&el);
        }

        let x_challenge_index_h1 = Verifier::<F, T, I, P, PR, FRI, PerRegisterARP>::bytes_to_challenge_index(
            &transcript.get_challenge_bytes(), 
            h1_lde.size(),
            self.lde_factor
        );

        let x_challenge_index_h2 = Verifier::<F, T, I, P, PR, FRI, PerRegisterARP>::bytes_to_challenge_index(
            &transcript.get_challenge_bytes(), 
            h2_lde.size(),
            self.lde_factor
        );


        let h1_fri_proof = FRI::prototype_into_proof(h1_fri_proof_proto, &h1_lde, x_challenge_index_h1)?;
        let h2_fri_proof = FRI::prototype_into_proof(h2_fri_proof_proto, &h2_lde, x_challenge_index_h2)?;


        let mut f_queries = vec![];
        for (o, lde) in f_oracles.iter().zip(f_ldes.iter()) {
            f_queries.push(o.query(x_challenge_index_h1, lde.as_ref()));
        }

        let g_query = g_oracle.query(x_challenge_index_h2, g_lde.as_ref());

        let proof = InstanceProof::<F, T, I, P, PR, FRI, PerRegisterARP> {
            f_at_z_m: f_at_z_m, 
            f_iop_roots: f_iop_roots,
            g_iop_root: g_iop_root,

            f_queries: f_queries,
            g_query: g_query,

            h1_iop_roots: h1_iop_roots,
            h2_iop_roots: h2_iop_roots,

            fri_proof_h1: h1_fri_proof,
            fri_proof_h2: h2_fri_proof,

            _marker_a: std::marker::PhantomData,
            _marker_t: std::marker::PhantomData,
            _marker_p: std::marker::PhantomData,
            _marker_fri: std::marker::PhantomData,
        };

        Ok(proof)
    }
}

#[test]
fn test_fib_prover() {
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

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();

    let lde_factor = 16;
    let output_at_degree_plus_one = 1;

    let prover = Prover::<Fr, Blake2sTranscript<Fr>, TrivialBlake2sIOP<Fr>, FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>, FRIProof<Fr, TrivialBlake2sIOP<Fr>>, NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>, PerRegisterARP>::new(
        props.clone(), 
        lde_factor,
        output_at_degree_plus_one
    ).expect("must work");

    let witness = witness.expect("some witness");

    let proof = prover.prove(witness).expect("must work");

    let verifier = Verifier::<Fr, Blake2sTranscript<Fr>, TrivialBlake2sIOP<Fr>, FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>, FRIProof<Fr, TrivialBlake2sIOP<Fr>>, NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>, PerRegisterARP>::new(
        props, 
        lde_factor
    ).expect("some verifier");

    println!("Verifier starts");
    let valid = verifier.verify(&proof).expect("must work");

    assert!(valid);
}

#[test]
fn test_soundness_of_fib_prover() {
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

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();

    let mut witness = witness.expect("some witness");

    witness[0][1] = Fr::from_str("123").unwrap();

    let lde_factor = 16;
    let output_at_degree_plus_one = 1;

    let prover = Prover::<Fr, Blake2sTranscript<Fr>, TrivialBlake2sIOP<Fr>, FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>, FRIProof<Fr, TrivialBlake2sIOP<Fr>>, NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>, PerRegisterARP>::new(
        props.clone(), 
        lde_factor,
        output_at_degree_plus_one
    ).expect("must work");

    let proof = prover.prove(witness).expect("must work");

    let verifier = Verifier::<Fr, Blake2sTranscript<Fr>, TrivialBlake2sIOP<Fr>, FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>, FRIProof<Fr, TrivialBlake2sIOP<Fr>>, NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>, PerRegisterARP>::new(
        props, 
        lde_factor
    ).expect("some verifier");

    println!("Verifier starts");
    let valid = verifier.verify(&proof).expect("must work");

    assert!(!valid);
}