use ff::*;

use crate::air::*;
use crate::fft::multicore::Worker;
use crate::ali::*;
use crate::arp::*;
use crate::iop::*;
use crate::iop::blake2s_trivial_iop::Blake2sIopTree;

use super::Fr;

pub struct DAStark<F: PrimeField> {
    pub n: usize,
    pub half_of_the_elements: Vec<F>,
}

impl<F: PrimeField> IntoARP<F> for DAStark<F> {
    fn into_arp(self) -> (Option<Vec<Vec<F>>>, InstanceProperties<F>) {
        assert!(self.half_of_the_elements.len() * 2 == self.n);

        let trace_length = 2 * 128 * self.n;

        let mut k = vec![];
        for i in 1..=127 {
            let el = F::from_str(&i.to_string()).unwrap();
            k.push(el);
        }

        // q^2 = B

        // (x + qy)^3 = x^3 + (q^2)q*y^3 + 3*x^2*q*y + 3x*(q^2)*y^2 = 
        // (x^3 + 3x*y^2*B + q(y^3*B + 3x^2*y*B))


        fn mimc_round<F: PrimeField>(el: (F, F), non_residue: &F, k: &[F], three: &F) -> Vec<(F, F)> {
            let mut trace = vec![];
            trace.push(el);

            let (mut c0, mut c1) = el;
            for k in k.iter() {
                let mut c0_squared = c0;
                c0_squared.square();

                let mut c0_cubed = c0_squared;
                c0_cubed.mul_assign(&c0);

                let mut c1_squared = c1;
                c1_squared.square();

                let mut c1_cubed = c1_squared;
                c1_cubed.mul_assign(&c1);

                let mut tmp_x = c0;
                tmp_x.mul_assign(&three);
                tmp_x.mul_assign(&non_residue);
                tmp_x.mul_assign(&c1_squared);

                let mut tmp_y = c1;
                tmp_x.mul_assign(&three);
                tmp_x.mul_assign(&non_residue);
                tmp_x.mul_assign(&c0_squared);

                let mut new_c0 = c0_cubed;
                new_c0.add_assign(&tmp_x);
                new_c0.add_assign(&k);

                let mut new_c1 = c1_cubed;
                new_c1.add_assign(&tmp_y);

                trace.push((new_c0, new_c1));

                c0 = new_c0;
                c1 = new_c1;
            }

            trace
        }

        fn mimc_hash<F: PrimeField>(xy: (F, F), non_residue: &F, k: &[F], three: &F) -> Vec<(F, F)> {
            let (x, y) = xy;
            
            let mut trace = vec![];
            let into_round = (x, F::zero());
            
            let round_0 = mimc_round(into_round, non_residue, k, three);

            let (_, r) = *round_0.last().unwrap();

            let into_round = (y, r);

            let round_1 = mimc_round(into_round, non_residue, k, three);

            trace.extend(round_0);
            trace.extend(round_1);

            trace
        }

        let mut non_residue = F::one();
        non_residue.negate();

        let num_registers = 4;

        let pc_register = Register::Register(0);
        let x_register = Register::Register(1);
        let y_register = Register::Register(2);
        let k_register = Register::Register(3);

        let pc_register_now = UnivariateTerm::<F>::from(pc_register);
        let x_register_now = UnivariateTerm::<F>::from(x_register);
        let y_register_now = UnivariateTerm::<F>::from(y_register);
        let k_register_now = UnivariateTerm::<F>::from(k_register);

        let mut pc_register_next_step = UnivariateTerm::<F>::from(pc_register);
        pc_register_next_step.set_step_difference(1);

        let mut x_register_next_step = UnivariateTerm::<F>::from(x_register);
        x_register_next_step.set_step_difference(1);

        let mut y_register_next_step = UnivariateTerm::<F>::from(y_register);
        y_register_next_step.set_step_difference(1);

        let mut x_register_next_step_squared = x_register_next_step.pow(2u64);
        let mut x_register_next_step_cubed = x_register_next_step.pow(3u64);

        let mut y_register_next_step_squared = y_register_next_step.pow(2u64);
        let mut y_register_next_step_cubed = y_register_next_step.pow(3u64);
        y_register_next_step_cubed *= non_residue;

        let three = F::from_str("3").unwrap();

        let mut x_cross_term = PolyvariateTerm::<F>::default();
        x_cross_term *= &three;
        x_cross_term *= x_register_next_step;
        x_cross_term *= y_register_next_step_squared;
        x_cross_term *= non_residue;

        let mut y_cross_term = PolyvariateTerm::<F>::default();
        y_cross_term *= &three;
        y_cross_term *= y_register_next_step;
        y_cross_term *= x_register_next_step_squared;
        y_cross_term *= non_residue;

        let mut constraint_x_mimc_round = Constraint::default();
        constraint_x_mimc_round.density = Box::from(DenseConstraint::default());

        constraint_x_mimc_round += x_register_next_step_cubed;
        constraint_x_mimc_round += x_cross_term;
        constraint_x_mimc_round += k_register_now;
        constraint_x_mimc_round -= x_register_now;

        let mut constraint_y_mimc_round = Constraint::default();
        constraint_y_mimc_round.density = Box::from(DenseConstraint::default());

        constraint_y_mimc_round += y_register_next_step_cubed;
        constraint_y_mimc_round += y_cross_term;
        constraint_y_mimc_round -= y_register_now;

        let mut zero_y_constraint = Constraint::default();
        zero_y_constraint.density = Box::from(SubdomainConstraint {
            subdomain_factor: 256,
            subdomain_coset_number: 0,
        });

        zero_y_constraint += y_register_now;

        let mut copy_y_constraint = Constraint::default();
        copy_y_constraint.density = Box::from(SubdomainConstraint {
            subdomain_factor: 256,
            subdomain_coset_number: 128,
        });

        copy_y_constraint += y_register_now;
        copy_y_constraint -= y_register_next_step;

        // while step is less than 128*n we would also use copy constraint

        let end_of_copy_constraint = 128*self.n;

        let mut pc_register_copy_constraint = Constraint::default();
        pc_register_copy_constraint.density = Box::from(AnywhereButInSetConstraint::new(vec![end_of_copy_constraint - 1]));
        pc_register_copy_constraint += pc_register_now;
        pc_register_copy_constraint -= pc_register_next_step;

        let mut pc_register_bump_constraint = Constraint::default();
        pc_register_bump_constraint.density = Box::from(DiscreteSetConstraint::new(vec![end_of_copy_constraint - 1]));
        pc_register_bump_constraint += pc_register_now;
        pc_register_bump_constraint -= pc_register_next_step;
        pc_register_bump_constraint -= &F::one();

        // we may also need one boundary constraint here, but skip for now

        let mut copy_x_constraint = Constraint::<F>::default();
        copy_x_constraint.density = Box::from(DenseConstraint::default());

        // hm, one can not easily restrict x(i) = x(2*i - 127) cause masking constraints can only have
        // additive shift
        let mut x_register_in_128 = UnivariateTerm::<F>::from(x_register);
        x_register_in_128.set_step_difference(127);

        // 😞😞😞😞😞😞😞😞😞😞😞😞😞😞😞😞😞

        unimplemented!()
    }
}