use ff::PrimeField;
use crate::polynomials::*;
use crate::arp::WitnessPolynomial;
use crate::arp::ARP;
use crate::air::*;
use crate::fft::multicore::Worker;

use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

impl<F: PrimeField> Hash for StepDifference<F> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match &self {
            StepDifference::Steps(s) => {
                s.hash(state);
            },
            StepDifference::Mask(m) => {
                m.into_raw_repr().as_ref().hash(state);
            }
        }
    }
}

// ARP works with remapped registers and no longer cares about their meaning
pub struct ALI<F: PrimeField> {
    pub f_poly: WitnessPolynomial<F>,
    pub g_poly: Polynomial<F, Values>,
}

impl<F: PrimeField> From<ARP<F>> for ALI<F> {
    fn from(arp: ARP<F>) -> ALI<F> {
        let mut max_constraint_power: u64 = 0;
        for c in arp.constraints.iter() {
            if c.degree > max_constraint_power {
                max_constraint_power = c.degree;
            }
        }

        // perform masking substitutions first 
        let mut all_masks = HashSet::<StepDifference<F>, _>::new();
        let mut mask_applied_polynomials = HashMap::<StepDifference<F>, Polynomial<F, Coefficients>, _>::new();

        fn get_masks_from_constraint<F: PrimeField>(
            set: &mut HashSet<StepDifference<F>>,
            constraint: &Constraint<F>
        ) {
            for t in constraint.terms.iter() {
                get_masks_from_term(set, t);
            }
        }

        fn get_masks_from_term<F: PrimeField>(
            set: &mut HashSet<StepDifference<F>>,
            constraint: &ConstraintTerm<F>
        ) {
            match constraint {
                ConstraintTerm::Univariate(uni) => {
                    set.insert(uni.steps_difference);
                    // if let StepDifference::Mask(m) = uni.steps_difference {
                    //     set.insert(m);
                    // }
                },
                ConstraintTerm::Polyvariate(poly) => {
                    for t in poly.terms.iter() {
                        set.insert(t.steps_difference);
                    }
                },
            }
        }

        for c in arp.constraints.iter() {
            get_masks_from_constraint(&mut all_masks, c);
        }

        fn evaluate_for_mask<F: PrimeField> (
            mut f: Polynomial<F, Coefficients>,
            mask: StepDifference<F>,
            worker: &Worker
        ) -> Polynomial<F, Coefficients> {
            match mask {
                StepDifference::Mask(m) => {
                    f.distribute_powers(&worker, m);
                },
                _ => {
                    unreachable!();
                }
            }

            f
        }

        let f = arp.witness_poly.expect("should me something");

        let worker = Worker::new();

        let f_poly = match &f {
            WitnessPolynomial::Single(p) => {
                p.clone()
            },
            _ => {
                unimplemented!();
            }
        };

        for mask in all_masks.into_iter() {
            mask_applied_polynomials.insert(mask, evaluate_for_mask(f_poly.clone(), mask, &worker));
        }

        let num_registers_sup = arp.num_registers.next_power_of_two();
        let num_steps_sup = arp.num_steps.next_power_of_two();

        assert_eq!(f_poly.as_ref().len(), num_registers_sup * num_steps_sup);

        let g_size = num_registers_sup * num_steps_sup * (max_constraint_power as usize);

        let g_poly = Polynomial::<F, Values>::new_for_size(g_size).expect("should work");

        ALI::<F> {
            f_poly: f,
            g_poly: g_poly
        }
    }
}

#[test]
fn test_fib_conversion() {
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::arp::IntoARP;
    use crate::ali::ALI;

    let fib = Fibonacci::<Fr> {
        first_a: Some(1),
        first_b: Some(1),
        final_a: Some(3),
        at_step: Some(2),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 5);
    let mut arp = ARP::<Fr>::new(test_tracer);
    arp.route_into_single_witness_poly().expect("must work");

    let ali = ALI::from(arp);
}