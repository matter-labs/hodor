use ff::PrimeField;
use crate::polynomials::*;
use crate::arp::WitnessPolynomial;
use crate::arp::ARP;
use crate::air::*;
use crate::fft::multicore::Worker;
use crate::SynthesisError;
use crate::domains::*;

pub mod deep_ali;

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
    pub mask_applied_polynomials: HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
    pub column_domain: Domain::<F>,
    pub full_trace_domain: Domain::<F>
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

        let boundary_constraint_mask = StepDifference::Mask(F::one());

        if all_masks.get(&boundary_constraint_mask).is_none() {
            all_masks.insert(boundary_constraint_mask);
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

        let column_domain = Domain::<F>::new_for_size(num_steps_sup as u64).expect("should be able to create");
        let full_trace_domain = Domain::<F>::new_for_size((num_steps_sup * num_registers_sup) as u64).expect("should be able to create");

        assert_eq!(f_poly.as_ref().len(), num_registers_sup * num_steps_sup);

        ALI::<F> {
            f_poly: f,
            g_poly: None,
            num_steps: arp.num_steps,
            num_registers: arp.num_registers,
            constraints: arp.constraints,
            max_constraint_power: max_constraint_power as usize,
            boundary_constraints: arp.boundary_constraints,
            mask_applied_polynomials: mask_applied_polynomials,
            column_domain,
            full_trace_domain
        }
    }
}

impl<F: PrimeField> ALI<F> {
    pub fn calculate_g(&mut self, alpha: F) -> Result<(), SynthesisError> {
        let mut current_coeff = F::one();
        // ---------------------

        // fn evaluate_constraint_term_into_values<F: PrimeField>(
        //     term: &ConstraintTerm<F>,
        //     substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
        //     worker: &Worker
        // ) -> Result<Polynomial<F, Values>, SynthesisError>
        // {
        //     let result = match term {
        //         ConstraintTerm::Univariate(uni) => {
        //             let mut base = substituted_witness.get(&uni.steps_difference).expect("should exist").clone();
        //             let scaling_factor = uni.power.next_power_of_two() as usize;
        //             base.extend(scaling_factor, &worker)?;
        //             let mut one = F::one();
   
        //             if uni.coeff != one {
        //                 one.negate();
        //                 if uni.coeff == one {
        //                     base.negate(&worker);
        //                 } else {
        //                     base.scale(&worker, uni.coeff);
        //                 }
        //             }

        //             let base = base.fft(&worker);

        //             base
        //         },
        //         ConstraintTerm::Polyvariate(_poly) => {
        //             unimplemented!();
        //         }
        //     };

        //     Ok(result)
        // }

        // ---------------------

        fn evaluate_constraint_term_into_coefficients<F: PrimeField>(
            term: &ConstraintTerm<F>,
            substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
            worker: &Worker
        ) -> Result<Polynomial<F, Coefficients>, SynthesisError>
        {
            let result = match term {
                ConstraintTerm::Univariate(uni) => {
                    let mut base = substituted_witness.get(&uni.steps_difference).expect("should exist").clone();
                    let mut new_base = if uni.power == 1u64 {
                        base
                    } else {
                        let scaling_factor = uni.power.next_power_of_two() as usize;
                        base.extend(scaling_factor, &worker)?;
                        let mut into_values = base.fft(&worker);
                        into_values.pow(&worker, uni.power);

                        into_values.ifft(&worker)
                    };

                    let one = F::one();
                    if uni.coeff != one {
                        let mut minus_one = one;
                        minus_one.negate();
                        if uni.coeff == minus_one {
                            new_base.negate(&worker);
                        } else {
                            new_base.scale(&worker, uni.coeff);
                        }
                    }

                    new_base
                },
                ConstraintTerm::Polyvariate(_poly) => {
                    unimplemented!();
                }
            };

            Ok(result)
        }

        // ---------------------

        let worker = Worker::new();
        let num_registers_sup = self.num_registers.next_power_of_two();
        let num_steps_sup = self.num_steps.next_power_of_two();
        let g_size = num_registers_sup * num_steps_sup * self.max_constraint_power;

        let mut g_poly = Polynomial::<F, Coefficients>::new_for_size(g_size).expect("should work");
        let subterm_coefficients = Polynomial::<F, Coefficients>::new_for_size(g_size).expect("should work");

        // ---------------------

        // such calls most likely will have start at 0 and num_steps = domain_size - 1
        fn inverse_divisor_for_dense_constraint_in_coset<F: PrimeField> (
            column_domain: &Domain<F>,
            term_evaluation_domain: &Domain<F>,
            alpha: F,
            start_at: u64,
            num_steps: u64,
            worker: &Worker
        ) -> Result<(Polynomial<F, Values>, usize), SynthesisError> {
            let mut divisor_degree = column_domain.size as usize;
            let divisor_domain_size = column_domain.size;
            divisor_degree -= start_at as usize;
            divisor_degree -= (divisor_domain_size - num_steps) as usize;

            let roots = {
                let roots_generator = column_domain.generator;

                let mut roots = vec![];
                let mut root = F::one();
                for _ in 0..start_at {
                    roots.push(root);
                    root.mul_assign(&roots_generator);                
                }

                let mut root = roots_generator.pow([num_steps]);
                for _ in num_steps..divisor_domain_size {
                    roots.push(root);
                    root.mul_assign(&roots_generator);
                }

                roots
            };

            let roots_iter = roots.iter();

            let evaluation_domain_generator = term_evaluation_domain.generator;
            let multiplicative_generator = F::multiplicative_generator();

            // these are values at the coset
            let mut inverse_divisors = Polynomial::<F, Values>::new_for_size(term_evaluation_domain.size as usize)?;

            // prepare for batch inversion
            worker.scope(inverse_divisors.as_ref().len(), |scope, chunk| {
                for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
                    scope.spawn(move |_| {
                        let mut x = evaluation_domain_generator.pow([(i*chunk) as u64]);
                        x.mul_assign(&multiplicative_generator);
                        for v in inv_divis.iter_mut() {
                            *v = x.pow([divisor_domain_size]);
                            v.sub_assign(&F::one());
                        }
                    });
                }
            });

            // now polynomial is filled with X^T - 1, and need to be inversed

            inverse_divisors.batch_inversion(&worker);

            // now do the evaluation

            worker.scope(inverse_divisors.as_ref().len(), |scope, chunk| {
                for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
                    let roots_iter_outer = roots_iter.clone();
                    scope.spawn(move |_| {
                        let mut x = evaluation_domain_generator.pow([(i*chunk) as u64]);
                        x.mul_assign(&multiplicative_generator);
                        for v in inv_divis.iter_mut() {
                            let mut c_inverse_by_alpha = *v;
                            c_inverse_by_alpha.mul_assign(&alpha);
                            let mut d = c_inverse_by_alpha;
                            for root in roots_iter_outer.clone() {
                                // (X - root)
                                let mut tmp = x;
                                tmp.sub_assign(&root);
                                d.mul_assign(&tmp);
                            } 
                            // alpha / ( (X^T-1) / (X - 1)(X - omega)(...) ) = alpha * (X - 1)(X - omega)(...) / (X^T-1)
                            *v = d;

                            x.mul_assign(&evaluation_domain_generator);
                        }
                    });
                }
            });

            // // TODO: optimize into batch inversion

            // worker.scope(inverse_divisors.as_ref().len(), |scope, chunk| {
            //     for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
            //         let roots_iter_outer = roots_iter.clone();
            //         scope.spawn(move |_| {
            //             let mut x = evaluation_domain_generator.pow([(i*chunk) as u64]);
            //             x.mul_assign(&multiplicative_generator);
            //             for v in inv_divis.iter_mut() {
            //                 let mut c = x.pow([divisor_domain_size]);
            //                 c.sub_assign(&F::one());
            //                 let mut c_inverse_by_alpha = c.inverse().expect("must exist");
            //                 c_inverse_by_alpha.mul_assign(&alpha);
            //                 let mut d = c_inverse_by_alpha;
            //                 for root in roots_iter_outer.clone() {
            //                     // (X - root)
            //                     let mut tmp = x;
            //                     tmp.sub_assign(&root);
            //                     d.mul_assign(&tmp);
            //                 } 
            //                 // alpha / ( (X^T-1) / (X - 1)(X - omega)(...) ) = alpha * (X - 1)(X - omega)(...) / (X^T-1)
            //                 *v = d;

            //                 x.mul_assign(&evaluation_domain_generator);
            //             }
            //         });
            //     }
            // });

            Ok((inverse_divisors, divisor_degree))
        }

        // ---------------------

        let subterm_domain = Domain::new_for_size(subterm_coefficients.as_ref().len() as u64)?;

        for constraint in self.constraints.iter() {
            current_coeff.mul_assign(&alpha);
            
            let mut subterm_coefficients = subterm_coefficients.clone();

            // first we need to calculate denominator at the value
            let (inverse_divisors, _divisor_degree) = match constraint.density {
                ConstraintDensity::Dense => {
                    let start_at = constraint.start_at as u64;

                    let result = inverse_divisor_for_dense_constraint_in_coset(
                        &self.column_domain,
                        &subterm_domain, 
                        current_coeff,
                        start_at,
                        self.num_steps as u64,
                        &worker
                    )?;

                    // println!("Divisor = {:?}", result);
                    result
                },
                _ => {
                    unimplemented!();
                }
            };

            // println!("Inverse divisors = {:?}", inverse_divisors);

            // println!("Evaluating constraint {:?}", constraint);

            for term in constraint.terms.iter() {
                let mut evaluated_term = evaluate_constraint_term_into_coefficients(
                    &term, 
                    &self.mask_applied_polynomials,
                    &worker
                )?;
                let factor = subterm_coefficients.as_ref().len() / evaluated_term.as_ref().len();
                evaluated_term.extend(factor, &worker)?;
                subterm_coefficients.add_assign(&worker, &evaluated_term);
            }

            subterm_coefficients.as_mut()[0].add_assign(&constraint.constant_term);

            // these values are correct and are evaluations of some polynomial at points (gen, gen * omega, gen * omega*2)
            let mut subterm_values_in_coset = subterm_coefficients.coset_fft(&worker);

            subterm_values_in_coset.mul_assign(&worker, &inverse_divisors);

            let subterm_coefficients = subterm_values_in_coset.icoset_fft(&worker);

            // let mut degree = subterm_coefficients.as_ref().len() - 1;
            // for c in subterm_coefficients.as_ref().iter().rev() {
            //     if c.is_zero() {
            //         degree -= 1;
            //     } else {
            //         break;
            //     }
            // }
            // println!("Final degree = {}", degree);

            g_poly.add_assign(&worker, &subterm_coefficients);
        }

        

        fn evaluate_boundary_constraint<F: PrimeField>(
            b_constraint: &BoundaryConstraint<F>,
            substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>
        ) -> Result<Polynomial<F, Coefficients>, SynthesisError>
        {
            let boundary_constraint_mask = StepDifference::Mask(F::one());
            let mut result = substituted_witness.get(&boundary_constraint_mask).expect("is some").clone();
            result.as_mut()[0].sub_assign(&b_constraint.value.expect("is some"));

            Ok(result)
        }

        for b_constraint in self.boundary_constraints.iter() {
            current_coeff.mul_assign(&alpha);

            // x - a
            let column_generator = self.column_domain.generator;
            let trace_generator = self.full_trace_domain.generator;

            let full_trace_size = self.full_trace_domain.size;
            let mut q_poly = Polynomial::<F, Coefficients>::new_for_size(full_trace_size as usize)?;
            q_poly.as_mut()[1] = F::one();
            let mut root = column_generator.pow([b_constraint.at_step as u64]);
            let reg_num = match b_constraint.register {
                Register::Register(reg_number) => {
                    reg_number
                },
                _ => {
                    unreachable!();
                }
            };

            root.mul_assign(&trace_generator.pow([reg_num as u64]));
            // omega^(t*W + i)
            q_poly.as_mut()[0].sub_assign(&root);

            let mut subterm_coefficients = subterm_coefficients.clone();

            let mut evaluated_term = evaluate_boundary_constraint(
                &b_constraint, 
                &self.mask_applied_polynomials
            )?;

            let factor = subterm_coefficients.as_ref().len() / evaluated_term.as_ref().len();
            evaluated_term.extend(factor, &worker)?;
            subterm_coefficients.add_assign(&worker, &evaluated_term);

            // TODO: optimize extensions for 1st degree polynomial

            let factor = subterm_coefficients.as_ref().len() / q_poly.as_ref().len();
            q_poly.extend(factor, &worker)?;

            let mut inverse_q_poly_coset_values = q_poly.coset_fft(&worker);
            inverse_q_poly_coset_values.batch_inversion(&worker);
            inverse_q_poly_coset_values.scale(&worker, current_coeff);

            // now those are in a form alpha * Q^-1

            let mut subterm_values_in_coset = subterm_coefficients.coset_fft(&worker);

            subterm_values_in_coset.mul_assign(&worker, &inverse_q_poly_coset_values);

            let subterm_coefficients = subterm_values_in_coset.icoset_fft(&worker);

            g_poly.add_assign(&worker, &subterm_coefficients);
        }


        self.g_poly = Some(g_poly);

        Ok(())
    }
}

#[test]
fn test_fib_conversion_into_ali() {
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::arp::IntoARP;
    use crate::ali::ALI;
    use crate::fft::multicore::Worker;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let mut arp = ARP::<Fr>::new(test_tracer);
    arp.route_into_single_witness_poly().expect("must work");

    let mut ali = ALI::from(arp);
    // println!("Mask applied polys = {:?}", ali.mask_applied_polynomials);
    let alpha = Fr::from_str("123").unwrap();
    ali.calculate_g(alpha).expect("must work");

    let g_poly_interpolant = ali.g_poly.take().expect("is something");
    println!("G coefficients = {:?}", g_poly_interpolant);
    let worker = Worker::new();
    let g_values = g_poly_interpolant.fft(&worker);
    println!("G values = {:?}", g_values);
}