use ff::PrimeField;
use crate::polynomials::*;
use crate::arp::WitnessPolynomial;
use crate::arp::ARP;
use crate::air::*;
use crate::fft::multicore::Worker;
use crate::SynthesisError;
use crate::domains::*;

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
#[derive(Debug)]
pub struct ALI<F: PrimeField> {
    pub f_poly: WitnessPolynomial<F>,
    pub g_poly: Option<Polynomial<F, Coefficients>>,
    pub num_steps: usize,
    pub num_registers: usize,
    pub max_constraint_power: usize,
    pub constraints: Vec<Constraint<F>>,
    pub boundary_constraints: Vec<BoundaryConstraint<F>>,
    mask_applied_polynomials: HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>
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

        ALI::<F> {
            f_poly: f,
            g_poly: None,
            num_steps: arp.num_steps,
            num_registers: arp.num_registers,
            constraints: arp.constraints,
            max_constraint_power: max_constraint_power as usize,
            boundary_constraints: arp.boundary_constraints,
            mask_applied_polynomials: mask_applied_polynomials
        }
    }
}

impl<F: PrimeField> ALI<F> {
    pub fn calculate_g(&mut self, alpha: F) -> Result<(), SynthesisError> {
        let mut current_coeff = F::one();
        let generator = F::multiplicative_generator();
        let column_domain = Domain::<F>::new_for_size(self.num_steps as u64)?;

        fn evaluate_constraint_term_in_coset<F: PrimeField>(
            term: &ConstraintTerm<F>,
            substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
            worker: &Worker
        ) -> Result<Polynomial<F, Values>, SynthesisError>
        {
            let result = match term {
                ConstraintTerm::Univariate(uni) => {
                    let mut base = substituted_witness.get(&uni.steps_difference).expect("should exist").clone();
                    let scaling_factor = uni.power.next_power_of_two() as usize;
                    base.extend(scaling_factor)?;
                    if uni.coeff != F::one() {
                        base.scale(&worker, uni.coeff);
                    }
                    let base = base.coset_fft(&worker);

                    base
                },
                ConstraintTerm::Polyvariate(_poly) => {
                    unimplemented!();
                }
            };

            Ok(result)
        }

        let worker = Worker::new();
        let num_registers_sup = self.num_registers.next_power_of_two();
        let num_steps_sup = self.num_steps.next_power_of_two();
        let g_size = num_registers_sup * num_steps_sup * self.max_constraint_power;

        let mut g_poly = Polynomial::<F, Values>::new_for_size(g_size).expect("should work");
        let subterm = g_poly.clone();

        // TODO: Check that witness values are evaluated at the coset, so division is valid

        for constraint in self.constraints.iter() {
            current_coeff.mul_assign(&alpha);
            // first we need to calculate denominator at the value
            let demonimator = match constraint.density {
                ConstraintDensity::Dense => {
                    let start_at = constraint.start_at;

                    let mut constraint_roots_generator = column_domain.generator;
                    constraint_roots_generator.mul_assign(&F::multiplicative_generator());

                    let mut root = constraint_roots_generator.pow([start_at as u64]);
                    let mut result = F::one();
                    for _ in start_at..self.num_steps {
                        let mut term = generator;
                        term.sub_assign(&root);
                        result.mul_assign(&term);
                        root.mul_assign(&constraint_roots_generator);
                    }

                    result
                },
                _ => {
                    unimplemented!();
                }
            };

            let demonimator = demonimator.inverse().expect("is non-zero");

            let mut subterm = subterm.clone();

            for term in constraint.terms.iter() {
                let mut evaluated_term = evaluate_constraint_term_in_coset(
                    &term, 
                    &self.mask_applied_polynomials,
                    &worker
                )?;

                let factor = subterm.as_ref().len() / evaluated_term.as_ref().len();

                evaluated_term.extend(factor)?;

                subterm.add_assign(&worker, &evaluated_term);
            }

            let mut alpha_by_demon = demonimator;
            alpha_by_demon.mul_assign(&current_coeff);

            // TODO: join two methods
            subterm.scale(&worker, alpha_by_demon);
            g_poly.add_assign(&worker, &subterm);
        }

        let g_interpolant = g_poly.icoset_fft(&worker);

        self.g_poly = Some(g_interpolant);

        Ok(())
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
    use crate::fft::multicore::Worker;

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

    let mut ali = ALI::from(arp);
    let alpha = Fr::from_str("3").unwrap();
    ali.calculate_g(alpha).expect("must work");

    let g_poly_interpolant = ali.g_poly.take().expect("is something");
    let worker = Worker::new();
    let g_values = g_poly_interpolant.fft(&worker);
    println!("G values = {:?}", g_values);
}