use ff::*;

use crate::air::*;
use crate::fft::multicore::Worker;
use crate::ali::*;
use crate::arp::*;
use crate::iop::*;
use crate::iop::blake2s_trivial_iop::Blake2sIopTree;

use super::Fr;

pub struct VDF<F: PrimeField> {
    pub start_c0: F,
    pub start_c1: F,
    pub num_operations: usize
}

impl<F: PrimeField> IntoARP<F> for VDF<F> {
    fn into_arp(self) -> (Option<Vec<Vec<F>>>, InstanceProperties<F>) {
        fn square<F: PrimeField>(el: (F, F), non_residue: &F) -> (F, F) {
            let (mut c0, mut c1) = el;
            let mut two_c0_c1 = c0;
            two_c0_c1.mul_assign(&c1);
            two_c0_c1.double();


            c0.square();
            c1.square();
            c1.mul_assign(&non_residue);
            c0.add_assign(&c1);

            (c0, two_c0_c1)
        }

        let mut non_residue = F::one();
        non_residue.negate();
        // TODO: check

        // squaring of (c0, c1) is (c0^2 + r*c1^2, 2*c0*c1)

        let num_registers = 2;

        let c0_register = Register::Register(0);
        let c1_register = Register::Register(1);

        let c0_value_now = UnivariateTerm::<F>::from(c0_register);
        let c1_value_now = UnivariateTerm::<F>::from(c1_register);

        let mut c0_value_next_step = UnivariateTerm::<F>::from(c0_register);
        c0_value_next_step.set_step_difference(1);

        let mut c1_value_next_step = UnivariateTerm::<F>::from(c1_register);
        c1_value_next_step.set_step_difference(1);

        let c0_squared = c0_value_now.pow(2u64);

        let mut c1_squared_by_r = c1_value_now.pow(2u64);
        c1_squared_by_r *= &non_residue;

        let mut two = F::one();
        two.double();

        let mut two_c0_c1 = PolyvariateTerm::<F>::default();
        two_c0_c1 *= &two;
        two_c0_c1 *= c0_value_now;
        two_c0_c1 *= c1_value_now;

        let mut constraint_0 = Constraint::default();
        constraint_0.density = ConstraintDensity::default();

        constraint_0 -= c0_squared;
        constraint_0 -= c1_squared_by_r;
        constraint_0 += c0_value_next_step;

        let mut constraint_1 = Constraint::default();
        constraint_1.density = ConstraintDensity::default();

        constraint_1 -= two_c0_c1;
        constraint_1 += c1_value_next_step;

        let num_values = self.num_operations + 1;

        let mut c0_witness = vec![F::zero(); num_values];
        c0_witness[0] = self.start_c0;
        let mut c1_witness = vec![F::zero(); num_values];
        c1_witness[0] = self.start_c1;

        let mut v0 = self.start_c0;
        let mut v1 = self.start_c1;
        for i in 0..self.num_operations {
            let tmp = square((v0, v1), &non_residue);
            v0 = tmp.0;
            v1 = tmp.1;
            c0_witness[i+1] = v0;
            c1_witness[i+1] = v1;
        }

        // add boundaty constraints
        let initial_c0_constraint = BoundaryConstraint::<F> {
            register: c0_register,
            at_row: 0,
            value: Some(self.start_c0)
        };
        let initial_c1_constraint = BoundaryConstraint::<F> {
            register: c1_register,
            at_row: 0,
            value: Some(self.start_c1)
        };

        let final_c0_constraint = BoundaryConstraint::<F> {
            register: c0_register,
            at_row: self.num_operations, // zero indexing
            value: c0_witness.last().cloned()
        };

        let final_c1_constraint = BoundaryConstraint::<F> {
            register: c1_register,
            at_row: self.num_operations, // zero indexing
            value: c1_witness.last().cloned()
        };

        let props = InstanceProperties::<F> {
            num_rows: self.num_operations + 1,
            num_registers: num_registers,
            constraints: vec![constraint_0, constraint_1],
            boundary_constraints: vec![initial_c0_constraint, initial_c1_constraint, final_c0_constraint, final_c1_constraint]
        };

        (Some(vec![c0_witness, c1_witness]), props)
    }
}

#[test]
fn try_prove_quadratic_vdf() {
    use std::time::Instant;
    use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;
    use crate::fri::*;
    use crate::transcript::*;
    use crate::ali::per_register::*;

    let vdf_instance = VDF::<Fr> {
        start_c0: Fr::one(),
        start_c1: Fr::one(),
        num_operations: (1 << 20) - 1
    };

    let worker = Worker::new();

    let lde_factor = 16;

    let mut transcript = Blake2sTranscript::new();

    let total_start = Instant::now();

    let start = Instant::now();
    let (witness, props) = vdf_instance.into_arp();
    println!("Done preraping and calculating VFD in {} ms", start.elapsed().as_millis());

    let witness = witness.expect("some witness");

    // let is_satisfied = ARPInstance::<Fr, PerRegisterARP>::is_satisfied(&props, &witness, &worker);
    // assert!(is_satisfied.is_ok());

    let arp = ARPInstance::<Fr, PerRegisterARP>::from_instance(props, &worker).expect("must work");

    let start = Instant::now();
    let witness_polys = arp.calculate_witness_polys(witness, &worker).expect("must work");
    println!("Witness polys taken {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let f_ldes: Vec<_> = witness_polys.iter().map(|w| {
        w.clone().lde(&worker, lde_factor).expect("must work")
    }).collect();
    println!("F LDEs is done after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let f_oracles: Vec<_> = f_ldes.iter().map(|l|
        Blake2sIopTree::create(l.as_ref())
    ).collect(); 
    println!("F oracles is done after {} ms", start.elapsed().as_millis());

    for o in f_oracles.iter() {
        transcript.commit_bytes(&o.get_root()[..]);
    }

    let start = Instant::now();
    let ali = ALIInstance::from_arp(arp, &worker).expect("is some");
    println!("ALI prepares after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let g_poly_interpolant = ali.calculate_g(&mut transcript, witness_polys.clone(), &worker).expect("is some");
    println!("G poly after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let g_lde = g_poly_interpolant.clone().lde(&worker, lde_factor).expect("is something");
    println!("G LDE is done after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let g_oracle = Blake2sIopTree::create(g_lde.as_ref());
    let z = g_oracle.get_challenge_scalar_from_root();
    println!("G oracle is done after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let (h1_lde, h2_lde) = ali.calculate_deep(
        &witness_polys,
        &f_ldes,
        &g_poly_interpolant,
        &g_lde,
        &mut transcript,
        &worker
    ).expect("must work");
    println!("H1 and H2 oracles done after {} ms", start.elapsed().as_millis());

    println!("Total proving time w/o FRI is {} ms", total_start.elapsed().as_millis());

    let h1_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h1_lde, lde_factor, 1, &worker);
    let h2_fri_proof = FRIIOP::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h2_lde, lde_factor, 1, &worker);

    println!("Total proving time with FRI is {} ms", total_start.elapsed().as_millis());
}