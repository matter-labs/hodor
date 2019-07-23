use ff::*;

use crate::air::*;
use crate::fft::multicore::Worker;
use crate::ali::*;
use crate::arp::*;
use crate::iop::*;
use crate::iop::blake2s_trivial_iop::Blake2sIopTree;

#[derive(PrimeField)]
#[PrimeFieldModulus = "3618502788666131213697322783095070105623107215331596699973092056135872020481"]
#[PrimeFieldGenerator = "7"]
pub struct Fr(FrRepr);

pub struct CubicVDF<F: PrimeField> {
    pub start_c0: F,
    pub start_c1: F,
    pub num_operations: usize
}

impl<F: PrimeField> IntoARP<F> for CubicVDF<F> {
    fn into_arp(self) -> ARP<F> {
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

        fn mul<F: PrimeField>(el: (F, F), by: &(F, F), non_residue: &F) -> (F, F) {
            let (c0, c1) = el;
            let (other_c0, other_c1) = by;

            let mut ac = c0;
            ac.mul_assign(other_c0);

            let mut bd_by_r = c1;
            bd_by_r.mul_assign(other_c1);
            bd_by_r.mul_assign(&non_residue);

            let mut new_c0 = ac;
            new_c0.add_assign(&bd_by_r);

            let mut ad = c0;
            ad.mul_assign(other_c1);

            let mut bc = c1;
            bc.mul_assign(other_c0);

            let mut new_c1 = ad;
            new_c1.add_assign(&bc);

            (new_c0, new_c1)
        }

        let mut non_residue = F::one();
        non_residue.negate();
        // TODO: check

        // squaring of (c0, c1) is (c0^2 + r*c1^2, 2*c0*c1)
        // multiplication is (a, b)*(c, d) = (a*c + r*b*d, a*c + b*d),
        // so we make four constraints of degree 2 and use 4 registers

        let num_registers = 4;

        let c0_register = Register::Register(0);
        let c1_register = Register::Register(1);
        let squaring_step_c0 = Register::Register(2);
        let squaring_step_c1 = Register::Register(3);

        let c0_value_now = UnivariateTerm::<F> {
            coeff: F::one(),
            register: c0_register,
            steps_difference: StepDifference::Steps(0),
            power: 1
        };

        let c1_value_now = UnivariateTerm::<F> {
            coeff: F::one(),
            register: c1_register,
            steps_difference: StepDifference::Steps(0),
            power: 1
        };

        let c0_value_next_step = UnivariateTerm::<F> {
            coeff: F::one(),
            register: c0_register,
            steps_difference: StepDifference::Steps(1),
            power: 1
        };

        let c1_value_next_step = UnivariateTerm::<F> {
            coeff: F::one(),
            register: c1_register,
            steps_difference: StepDifference::Steps(1),
            power: 1
        };

        let c0_squared = UnivariateTerm::<F> {
            coeff: F::one(),
            register: c0_register,
            steps_difference: StepDifference::Steps(0),
            power: 2
        };

        let c1_squared_by_r = UnivariateTerm::<F> {
            coeff: non_residue,
            register: c1_register,
            steps_difference: StepDifference::Steps(0),
            power: 2
        };

        let squaring_step_intermediate_c0 = UnivariateTerm::<F> {
            coeff: F::one(),
            register: squaring_step_c0,
            steps_difference: StepDifference::Steps(0),
            power: 1
        };

        let squaring_step_intermediate_c1 = UnivariateTerm::<F> {
            coeff: F::one(),
            register: squaring_step_c1,
            steps_difference: StepDifference::Steps(0),
            power: 1
        };

        let mut two = F::one();
        two.double();

        let two_c0_c1 = PolyvariateTerm::<F> {
            coeff: two,
            terms: vec![c0_value_now, c1_value_now],
            total_degree: 2u64
        };

        // intermediate c0 = c0^2 + r*c1^2 
        let mut constraint_0 = Constraint::default();
        constraint_0.start_at = 0;
        constraint_0.density = ConstraintDensity::Dense;

        constraint_0 -= c0_squared;
        constraint_0 -= c1_squared_by_r;
        constraint_0 += squaring_step_intermediate_c0;

        // same for intermediate c1

        let mut constraint_1 = Constraint::default();
        constraint_1.start_at = 0;
        constraint_1.density = ConstraintDensity::Dense;

        constraint_1 -= two_c0_c1;
        constraint_1 += squaring_step_intermediate_c1;

        // constraint new_c0 = c0*squaring_step_intermediate_c0 + r*c1*squaring_step_intermediate_c1

        let c0_by_intermediate = PolyvariateTerm::<F> {
            coeff: F::one(),
            terms: vec![c0_value_now.clone(), squaring_step_intermediate_c0.clone()],
            total_degree: 2u64
        };

        let r_c1_by_intermediate = PolyvariateTerm::<F> {
            coeff: non_residue,
            terms: vec![c1_value_now.clone(), squaring_step_intermediate_c1.clone()],
            total_degree: 2u64
        };

        let mut constraint_2 = Constraint::default();
        constraint_2.start_at = 0;
        constraint_2.density = ConstraintDensity::Dense;

        constraint_2 -= c0_by_intermediate;
        constraint_2 -= r_c1_by_intermediate;
        constraint_2 += c0_value_next_step;

        // constraint new_c1 = c0*squaring_step_intermediate_c1 + c1*squaring_step_intermediate_c0

        let c0_by_intermediate_c1 = PolyvariateTerm::<F> {
            coeff: F::one(),
            terms: vec![c0_value_now.clone(), squaring_step_intermediate_c1.clone()],
            total_degree: 2u64
        };

        let c1_by_intermediate_c0 = PolyvariateTerm::<F> {
            coeff: F::one(),
            terms: vec![c1_value_now.clone(), squaring_step_intermediate_c0.clone()],
            total_degree: 2u64
        };

        let mut constraint_3 = Constraint::default();
        constraint_3.start_at = 0;
        constraint_3.density = ConstraintDensity::Dense;

        constraint_3 -= c0_by_intermediate_c1;
        constraint_3 -= c1_by_intermediate_c0;
        constraint_3 += c1_value_next_step;

        let num_values = self.num_operations + 1;


        let mut c0_witness = vec![F::zero(); num_values];
        c0_witness[0] = self.start_c0;
        let mut c1_witness = vec![F::zero(); num_values];
        c1_witness[0] = self.start_c1;
        let mut c0_squaring_witness = vec![F::zero(); num_values];
        let mut c1_squaring_witness = vec![F::zero(); num_values];

        let tmp = square((self.start_c0, self.start_c1), &non_residue);

        c0_squaring_witness[0] = tmp.0;
        c1_squaring_witness[0] = tmp.1;


        let mut v0 = self.start_c0;
        let mut v1 = self.start_c1;
        let mut v0_squared = tmp.0;
        let mut v1_squared = tmp.1;
        for i in 0..self.num_operations {
            let cubed = mul((v0_squared, v1_squared), &(v0, v1), &non_residue);
            v0 = cubed.0;
            v1 = cubed.1;
            c0_witness[i+1] = v0;
            c1_witness[i+1] = v1;

            let tmp = square((v0, v1), &non_residue);
            v0_squared = tmp.0;
            v1_squared = tmp.1;
            c0_squaring_witness[i+1] = tmp.0;
            c1_squaring_witness[i+1] = tmp.1;
        }

        // add boundaty constraints
        let initial_c0_constraint = BoundaryConstraint::<F> {
            register: c0_register,
            at_step: 0,
            value: Some(self.start_c0)
        };
        let initial_c1_constraint = BoundaryConstraint::<F> {
            register: c1_register,
            at_step: 0,
            value: Some(self.start_c1)
        };

        let final_c0_constraint = BoundaryConstraint::<F> {
            register: c0_register,
            at_step: self.num_operations + 1,
            value: c0_witness.last().cloned()
        };

        let final_c1_constraint = BoundaryConstraint::<F> {
            register: c1_register,
            at_step: self.num_operations + 1,
            value: c1_witness.last().cloned()
        };

        ARP::<F> {
            witness: Some(vec![c0_witness, c1_witness, c0_squaring_witness, c1_squaring_witness]),
            witness_poly: None,
            num_steps: self.num_operations,
            num_registers: num_registers,
            constraints: vec![constraint_0, constraint_1, constraint_2, constraint_3],
            boundary_constraints: vec![initial_c0_constraint, initial_c1_constraint, final_c0_constraint, final_c1_constraint]
        }
    }
}

#[test]
fn try_prove_cubic_vdf() {
    use std::time::Instant;

    let vdf_instance = CubicVDF::<Fr> {
        start_c0: Fr::one(),
        start_c1: Fr::one(),
        num_operations: (1 << 20) - 1
    };

    let worker = Worker::new();

    let lde_factor = 16;

    let total_start = Instant::now();

    let start = Instant::now();

    let mut arp = ARP::new(vdf_instance);
    println!("Done preraping and calculating VFD in {} ms", start.elapsed().as_millis());

    let start = Instant::now();

    arp.route_into_single_witness_poly().expect("must work");

    let f_witness_interpolant = arp.witness_poly.clone().expect("is something");
    let f_witness_interpolant = match f_witness_interpolant {
        WitnessPolynomial::Single(f) => {
            f
        },
        _ => {
            unreachable!();
        }
    };

    let f_lde = f_witness_interpolant.lde(&worker, lde_factor).expect("is something");
    println!("F LDE is done after {} ms", start.elapsed().as_millis());
    let start = Instant::now();

    let f_oracle = Blake2sIopTree::create(f_lde.as_ref());
    println!("F oracle is done after {} ms", start.elapsed().as_millis());
    let start = Instant::now();

    let mut ali = ALI::from(arp);
    let alpha = f_oracle.get_challenge_scalar_from_root();
    // let alpha = Fr::from_str("123").unwrap();
    ali.calculate_g(alpha).expect("must work");

    println!("G is calculated after {} ms", start.elapsed().as_millis());
    let start = Instant::now();

    let g_poly_interpolant = ali.g_poly.clone().expect("is something");
    let g_lde = g_poly_interpolant.lde(&worker, lde_factor).expect("is something");

    println!("G LDE is done after {} ms", start.elapsed().as_millis());
    let start = Instant::now();

    let g_oracle = Blake2sIopTree::create(g_lde.as_ref());
    let z = g_oracle.get_challenge_scalar_from_root();
    // let z = Fr::from_str("62").unwrap();

    println!("G oracle is done after {} ms", start.elapsed().as_millis());
    let start = Instant::now();

    let mut deep_ali = DeepALI::from(ali);
    deep_ali.make_deep(f_lde, g_lde, z).expect("must work");

    println!("H1 and H2 oracles done after {} ms", start.elapsed().as_millis());
    println!("Total proving time w/o FRI is {} ms", total_start.elapsed().as_millis());
}