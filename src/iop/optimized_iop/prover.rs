use std::convert::TryInto;

use ff::{PrimeField, PrimeFieldRepr, Field};
use rand::{Rng, SeedableRng, XorShiftRng};

use crate::air::IntoAIR;
use crate::ali::per_register::*;
use crate::ali::*;
use crate::arp::*;
use crate::fft::multicore::Worker;
// use crate::fri::*;
use super::fri::*;
use super::iop::{IOP};
use super::query::IopQuery;
use crate::transcript::Transcript;
use crate::transcript::*;
use super::verifier::*;
use crate::SynthesisError;
use super::fri::{FriProofPrototype, FriProof};

/*

This module contains a stand-alone prover that is basically a composition from all other modules. Workflow is simple: new prover is constructed from the instance of the "problem"
and does all the precomputations for ARP and ALI step, e.g. divisor precomputation for ALI. It takes additional parameters
such as a desired LDE factor (1/rho paramter for FRI) and at what step FRI should output the final polynomial in a form
of plain coefficients

Test containts a use example

*/

pub struct Prover<
    F: PrimeField,
    T: Transcript<F>,
    I: IOP<F>,
    P: FriProofPrototype<F, I>,
    PR: FriProof<F, I>,
    FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR>,
    A: ARPType,
> {
    arp: ARPInstance<F, A>,
    ali: ALIInstance<F, A>,
    lde_factor: usize,
    fri_final_degree_plus_one: usize,
    worker: Worker,

    _marker_t: std::marker::PhantomData<T>,
    _marker_i: std::marker::PhantomData<I>,
    _marker_p: std::marker::PhantomData<P>,
    _marker_pr: std::marker::PhantomData<PR>,
    _marker_fri: std::marker::PhantomData<FRI>,
}

impl<
        F: PrimeField,
        T: Transcript<F>,
        I: IOP<F>,
        P: FriProofPrototype<F, I>,
        PR: FriProof<F, I>,
        FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR>,
    > Prover<F, T, I, P, PR, FRI, PerRegisterARP>
{
    pub fn new(
        instance: InstanceProperties<F>,
        lde_factor: usize,
        fri_final_degree_plus_one: usize,
    ) -> Result<Self, SynthesisError> {
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

    pub fn prove(
        &self,
        witness: Vec<Vec<F>>,
    ) -> Result<InstanceProof<F, T, I, P, PR, FRI, PerRegisterARP>, SynthesisError> {
        let mut transcript = T::new();

        let witness_polys = self.arp.calculate_witness_polys(witness, &self.worker)?;

        let mut f_ldes = Vec::with_capacity(witness_polys.len());

        for w in witness_polys.iter() {
            let f_lde = w.clone().lde(&self.worker, self.lde_factor)?;
            f_ldes.push(f_lde);
        }

        let number_of_ldes = f_ldes.len();
        let length_of_single_lde = f_ldes[0].as_ref().len();

        let mut combined_lde = vec![];
        for lde in f_ldes.iter(){
            combined_lde.extend_from_slice(lde.as_ref());
        }

        assert_eq!(combined_lde.len(), number_of_ldes*length_of_single_lde);

        let single_oracle_from_multiple_lde = I::create(&combined_lde, length_of_single_lde);
        // let q = single_oracle_from_multiple_lde.query(2, &combined_lde, length_of_single_lde);
        // assert_eq!(q.value().len(), number_of_ldes, "length of combined leaf does not match with number of ldes");
        let root_of_combined_tree = single_oracle_from_multiple_lde.get_root();
        transcript.commit_bytes(&root_of_combined_tree.as_ref());


        let g_poly_interpolant =
            self.ali
                .calculate_g(&mut transcript, witness_polys.clone(), &self.worker)?;

        let g_lde = g_poly_interpolant
            .clone()
            .lde(&self.worker, self.lde_factor)?;
        
        let g_oracle = I::create(g_lde.as_ref(), g_lde.as_ref().len());

        let g_iop_root = g_oracle.get_root();
        transcript.commit_bytes(g_iop_root.as_ref());

        

        // println!("Calculating DEEP part");

        let (h1_lde, h2_lde, f_at_z_m, _g_at_z) = self.ali.calculate_deep(
            &witness_polys,
            &f_ldes,
            &g_poly_interpolant,
            &g_lde,
            &mut transcript,
            &self.worker,
        )?;

        let fri_final_poly_degree = self.fri_final_degree_plus_one;

        // println!("Calculating FRI part");

        let h1_fri_proof_proto = FRI::proof_from_lde(
            &h1_lde,
            self.lde_factor,
            fri_final_poly_degree,
            &self.worker,
        )?;

        let h2_fri_proof_proto = FRI::proof_from_lde(
            &h2_lde,
            self.lde_factor,
            fri_final_poly_degree,
            &self.worker,
        )?;

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

        let x_challenge_index_h1 =
            Verifier::<F, T, I, P, PR, FRI, PerRegisterARP>::bytes_to_challenge_index(
                &transcript.get_challenge_bytes(),
                h1_lde.size(),
                self.lde_factor,
            );        

        let x_challenge_index_h2 =
            Verifier::<F, T, I, P, PR, FRI, PerRegisterARP>::bytes_to_challenge_index(
                &transcript.get_challenge_bytes(),
                h2_lde.size(),
                self.lde_factor,
            );            

        let h1_fri_proof =
            FRI::prototype_into_proof(h1_fri_proof_proto, &h1_lde, x_challenge_index_h1)?;

        let h2_fri_proof =
        FRI::prototype_into_proof(h2_fri_proof_proto, &h2_lde, x_challenge_index_h2)?;


        let combined_f_query = single_oracle_from_multiple_lde.query(x_challenge_index_h1, &combined_lde, length_of_single_lde);

        let g_query = g_oracle.query(x_challenge_index_h2, g_lde.as_ref(), g_lde.as_ref().len());

        let proof = InstanceProof::<F, T, I, P, PR, FRI, PerRegisterARP> {
            f_at_z_m: f_at_z_m,
            f_iop_root: root_of_combined_tree,
            g_iop_root: g_iop_root,

            f_query: combined_f_query,
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
fn test_optimized_fib_prover() {
    use crate::air::Fibonacci;
    use crate::air::IntoAIR;
    use crate::air::TestTraceSystem;
    use crate::ali::per_register::*;
    use crate::ali::*;
    use crate::arp::*;
    use crate::fft::multicore::Worker;
    // use crate::fri::*;
    use super::fri::*;
    // use crate::iop::blake2s_trivial_iop::Blake2sIopTree;
    // use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;
    use super::tree::Blake2sIopTree;
    // use super::TrivialBlake2sIOP;
    use super::iop::TrivialBlake2sIOP;
    use crate::transcript::Transcript;
    use crate::transcript::*;
    use crate::Fr;
    use ff::Field;
    use super::prover::Prover;
    use super::verifier::Verifier;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData,
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();

    let lde_factor = 16;
    let output_at_degree_plus_one = 1;

    let prover = Prover::<
        Fr,
        Blake2sTranscript<Fr>,
        TrivialBlake2sIOP<Fr>,
        FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>,
        FRIProof<Fr, TrivialBlake2sIOP<Fr>>,
        NaiveFriIop<Fr, TrivialBlake2sIOP<Fr>>,
        PerRegisterARP,
    >::new(props.clone(), lde_factor, output_at_degree_plus_one)
    .expect("must work");

    let witness = witness.expect("some witness");

    let proof = prover.prove(witness).expect("must work");

    let verifier = Verifier::<
        Fr,
        Blake2sTranscript<Fr>,
        TrivialBlake2sIOP<Fr>,
        FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>,
        FRIProof<Fr, TrivialBlake2sIOP<Fr>>,
        NaiveFriIop<Fr, TrivialBlake2sIOP<Fr>>,
        PerRegisterARP,
    >::new(props, lde_factor)
    .expect("some verifier");

    println!("Verifier starts");
    let valid = verifier.verify(&proof).expect("must work");

    assert!(valid);
}

#[test]
fn test_soundness_of_fib_prover() {
    use crate::air::Fibonacci;
    use crate::air::IntoAIR;
    use crate::air::TestTraceSystem;
    use crate::ali::per_register::*;
    use crate::ali::*;
    use crate::arp::*;
    use crate::fft::multicore::Worker;
    // use crate::fri::*;
    use super::fri::*;
    use super::tree::Blake2sIopTree;
    use super::iop::TrivialBlake2sIOP;
    // use crate::iop::blake2s_trivial_iop::Blake2sIopTree;
    // use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;
    use crate::transcript::Transcript;
    use crate::transcript::*;
    use crate::Fr;
    use ff::Field;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData,
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();

    let mut witness = witness.expect("some witness");

    witness[0][1] = Fr::from_str("123").unwrap();

    let lde_factor = 16;
    let output_at_degree_plus_one = 1;

    let prover = Prover::<
        Fr,
        Blake2sTranscript<Fr>,
        TrivialBlake2sIOP<Fr>,
        FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>,
        FRIProof<Fr, TrivialBlake2sIOP<Fr>>,
        NaiveFriIop<Fr, TrivialBlake2sIOP<Fr>>,
        PerRegisterARP,
    >::new(props.clone(), lde_factor, output_at_degree_plus_one)
    .expect("must work");

    let proof = prover.prove(witness).expect("must work");

    let verifier = Verifier::<
        Fr,
        Blake2sTranscript<Fr>,
        TrivialBlake2sIOP<Fr>,
        FRIProofPrototype<Fr, TrivialBlake2sIOP<Fr>>,
        FRIProof<Fr, TrivialBlake2sIOP<Fr>>,
        NaiveFriIop<Fr, TrivialBlake2sIOP<Fr>>,
        PerRegisterARP,
    >::new(props, lde_factor)
    .expect("some verifier");

    println!("Verifier starts");
    let valid = verifier.verify(&proof).expect("must work");

    assert!(!valid);
}


#[test]
fn test_multiple_fr_into_buf(){
    use crate::Fr;
    use crate::ff::{PrimeField, Field, PrimeFieldRepr};
    use rand::Rand;
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let mut values = vec![Fr::zero();5];

    for value in values.iter_mut(){
        *value = Fr::rand(&mut rng);
    }

    let mut buf = vec![0u8; values.len()*32];

    for (idx, el) in values.iter().enumerate(){
        let repr = el.into_raw_repr();
        repr.write_le(&mut buf[idx*32..(idx+1)*32][..]).expect("will write");
        println!("el {:?}", repr);
        println!("encoded {:?}", buf[idx*32..(idx+1)*32].to_vec());
    }

}