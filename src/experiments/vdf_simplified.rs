use ff::*;

use crate::air::*;
use crate::fft::multicore::Worker;
use crate::ali::*;
use crate::arp::*;
use crate::iop::*;
use crate::iop::blake2s_trivial_iop::Blake2sIopTree;

use super::fields::vdf_128_ef::Fr;

pub struct VDF<F: PrimeField> {
    pub x0: F,
    pub x1: F,
    pub num_operations: usize
}

impl<F: PrimeField> IntoARP<F> for VDF<F> {
    fn into_arp(self) -> (Option<Vec<Vec<F>>>, InstanceProperties<F>) {
        let num_rows = self.num_operations + 1;
        assert!(num_rows.is_power_of_two());
        // def round(x):
        //     return [cube_root(x[0] + x[1]), (x[1] + 1) % p]

        // so going backwards:

        // (x0 + x1)^3 = x0_next
        // x1 - 1 = x1_next

        // make it as 
        // (x0 + x1) ^2 = tmp
        // tmp * (x0 + x1) = x0_next
        // x1 - 1 = x1_next


        // x0^3 = x0 + x1
        // x1 = x1 - 1

        // (x0 + x1) ^3 = x0
        // 

        let num_registers = 3;

        let x0_register = Register::Register(0);
        let x1_register = Register::Register(1);
        let x0_plus_x1_squared_register = Register::Register(2);

        let x0_value_now = UnivariateTerm::<F>::from(x0_register);
        let x1_value_now = UnivariateTerm::<F>::from(x1_register);
        let tmp_value_now = UnivariateTerm::<F>::from(x0_plus_x1_squared_register);

        let mut x0_value_next_step = UnivariateTerm::<F>::from(x0_register);
        x0_value_next_step.set_step_difference(1);

        let mut x1_value_next_step = UnivariateTerm::<F>::from(x1_register);
        x1_value_next_step.set_step_difference(1);

        let x0_squared = x0_value_now.pow(2u64);
        let x1_squared = x1_value_now.pow(2u64);

        let one = F::one();

        let mut two = F::one();
        two.double();

        let mut two_x0_x1 = PolyvariateTerm::<F>::default();
        two_x0_x1 *= &two;
        two_x0_x1 *= x0_value_now;
        two_x0_x1 *= x1_value_now;

        let anywhere_but_last_density = AnywhereButInSetConstraint::new(vec![num_rows-1]);
        let anywhere_but_last_density = Box::from(anywhere_but_last_density);
        let mut constraint_0 = Constraint::default();
        // constraint_0.density = Box::from(DenseConstraint::default());
        constraint_0.density = anywhere_but_last_density.clone();
        constraint_0 += UnivariateTerm::<F>::from(tmp_value_now);
        constraint_0 -= x0_squared;
        constraint_0 -= x1_squared;
        constraint_0 -= two_x0_x1;

        let mut x0_by_tmp = PolyvariateTerm::<F>::default();
        x0_by_tmp *= x0_value_now;
        x0_by_tmp *= tmp_value_now;

        let mut x1_by_tmp = PolyvariateTerm::<F>::default();
        x1_by_tmp *= x1_value_now;
        x1_by_tmp *= tmp_value_now;

        let mut constraint_1 = Constraint::default();
        constraint_1.density = anywhere_but_last_density.clone();
        constraint_1 += x0_value_next_step;
        constraint_1 -= x0_by_tmp;
        constraint_1 -= x1_by_tmp;

        let mut constraint_2 = Constraint::default();
        constraint_2.density = anywhere_but_last_density.clone();
        constraint_2 += x1_value_next_step;
        constraint_2 -= x1_value_now;
        constraint_2 += one;

        let num_values = num_rows;

        let mut x0_witness = vec![F::zero(); num_values];
        x0_witness[0] = self.x0;
        let mut x1_witness = vec![F::zero(); num_values];
        x1_witness[0] = self.x1;
        let mut tmp_values = vec![F::zero(); num_values];

        let mut v0 = self.x0;
        let mut v1 = self.x1;
        for i in 0..self.num_operations {
            // square
            let mut tmp = v0;
            tmp.add_assign(&v1);
            tmp.square();
            tmp_values[i] = tmp;

            // cube
            let mut t0 = v0;
            t0.add_assign(&v1);
            t0.mul_assign(&tmp);

            let mut t1 = v1;
            t1.sub_assign(&F::one());

            v0 = t0;
            v1 = t1;
            x0_witness[i+1] = v0;
            x1_witness[i+1] = v1;
        }

        // add boundaty constraints
        let initial_x0_constraint = BoundaryConstraint::<F> {
            register: x0_register,
            at_row: 0,
            value: Some(self.x0)
        };

        let initial_x1_constraint = BoundaryConstraint::<F> {
            register: x1_register,
            at_row: 0,
            value: Some(self.x1)
        };

        let final_x0_constraint = BoundaryConstraint::<F> {
            register: x0_register,
            at_row: self.num_operations, // zero indexing
            value: x0_witness.last().cloned()
        };

        let final_x1_constraint = BoundaryConstraint::<F> {
            register: x1_register,
            at_row: self.num_operations, // zero indexing
            value: x1_witness.last().cloned()
        };

        let props = InstanceProperties::<F> {
            num_rows: self.num_operations + 1,
            num_registers: num_registers,
            constraints: vec![constraint_0, constraint_1, constraint_2],
            boundary_constraints: vec![initial_x0_constraint, initial_x1_constraint, final_x0_constraint, final_x1_constraint]
        };

        (Some(vec![x0_witness, x1_witness, tmp_values]), props)
    }
}

#[test]
fn try_prove_ef_vdf() {
    use std::time::Instant;
    use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;
    use crate::fri::*;
    use crate::transcript::*;
    use crate::ali::per_register::*;

    let mut two = Fr::one();
    two.double();
    let vdf_instance = VDF::<Fr> {
        x0: Fr::one(),
        x1: two,
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

    let is_satisfied = ARPInstance::<Fr, PerRegisterARP>::is_satisfied(&props, &witness, &worker);
    assert!(is_satisfied.is_ok(), "trace is invalid with error {}", is_satisfied.err().unwrap());

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
    let ali = ALIInstance::from_arp(arp, &worker).unwrap();
    println!("ALI prepares after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let g_poly_interpolant = ali.calculate_g(&mut transcript, witness_polys.clone(), &worker).expect("is some");
    println!("G poly after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let g_lde = g_poly_interpolant.clone().lde(&worker, lde_factor).expect("is something");
    println!("G LDE is done after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let g_oracle = Blake2sIopTree::create(g_lde.as_ref());
    transcript.commit_bytes(&g_oracle.get_root());
    println!("G oracle is done after {} ms", start.elapsed().as_millis());

    let start = Instant::now();
    let (h1_lde, h2_lde, _, _) = ali.calculate_deep(
        &witness_polys,
        &f_ldes,
        &g_poly_interpolant,
        &g_lde,
        &mut transcript,
        &worker
    ).expect("must work");
    println!("H1 and H2 oracles done after {} ms", start.elapsed().as_millis());

    println!("Total proving time w/o FRI is {} ms", total_start.elapsed().as_millis());

    let h1_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h1_lde, lde_factor, 1, &worker);
    let h2_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde_by_values(&h2_lde, lde_factor, 1, &worker);

    println!("Total proving time with FRI is {} ms", total_start.elapsed().as_millis());
}