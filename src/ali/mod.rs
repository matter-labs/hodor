use ff::PrimeField;
use crate::polynomials::*;
use crate::arp::WitnessPolynomial;
use crate::arp::ARP;
use crate::air::*;
use crate::fft::multicore::Worker;
use crate::SynthesisError;
use crate::domains::*;

pub mod deep_ali;

pub use self::deep_ali::*;

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
    pub inversed_divisors_in_cosets: HashMap::<ConstraintDensity, Polynomial<F, Values>>,
    pub column_domain: Domain::<F>,
    pub full_trace_domain: Domain::<F>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WitnessEvaluationData<F: PrimeField> {
    mask: StepDifference<F>,
    power: u64,
    total_lde_length: u64
}

impl<F: PrimeField> Hash for WitnessEvaluationData<F> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.mask.hash(state);
        self.power.hash(state);
        self.total_lde_length.hash(state);
    }
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

        assert_eq!(f_poly.size(), num_registers_sup * num_steps_sup);

        ALI::<F> {
            f_poly: f,
            g_poly: None,
            num_steps: arp.num_steps,
            num_registers: arp.num_registers,
            constraints: arp.constraints,
            max_constraint_power: max_constraint_power as usize,
            boundary_constraints: arp.boundary_constraints,
            mask_applied_polynomials: mask_applied_polynomials,
            inversed_divisors_in_cosets: HashMap::new(),
            column_domain,
            full_trace_domain
        }
    }
}

impl<F: PrimeField> ALI<F> {
    pub fn calculate_g(&mut self, alpha: F) -> Result<(), SynthesisError> {
        // ---------------------

        fn evaluate_constraint_term_into_coefficients<F: PrimeField>(
            term: &ConstraintTerm<F>,
            substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut HashMap::<WitnessEvaluationData<F>, Polynomial<F, Values>>,
            _power_hint: Option<u64>,
            worker: &Worker
        ) -> Result<Polynomial<F, Coefficients>, SynthesisError>
        {
            let result = match term {
                ConstraintTerm::Univariate(uni) => {
                    let t = evaluate_univariate_term_into_coefficients(
                        uni,
                        substituted_witness,
                        &worker 
                    )?;

                    t
                },
                ConstraintTerm::Polyvariate(poly) => {
                    let power_hint = Some(poly.total_degree as u64);
                    let mut values_result: Option<Polynomial<F, Values>> = None;
                    // evaluate subcomponents in a value form and multiply
                    for uni in poly.terms.iter() {
                        let t = evaluate_univariate_term_into_values(
                            uni, 
                            substituted_witness,
                            evaluated_univariate_terms,
                            power_hint,
                            &worker
                        )?;
                        if let Some(res) = values_result.as_mut() {
                            res.mul_assign(&worker, &t);
                        } else {
                            values_result = Some(t); 
                        }
                    }

                    let as_values = values_result.expect("is some");
                    // go back into coefficients form and scale
                    let mut as_coeffs = as_values.ifft(&worker);

                    as_coeffs.scale(&worker, poly.coeff);

                    as_coeffs
                }
            };

            Ok(result)
        }

        // ---------------------

        fn evaluate_univariate_term_into_coefficients<F: PrimeField>(
            uni: &UnivariateTerm<F>,
            substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
            worker: &Worker
        ) -> Result<Polynomial<F, Coefficients>, SynthesisError>
        {
            let base = substituted_witness.get(&uni.steps_difference).expect("should exist").clone();

            let mut new_base = if uni.power == 1u64 {
                base
            } else {
                let factor = uni.power.next_power_of_two() as usize;
                let mut into_values = base.lde(&worker, factor)?;
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

            Ok(new_base)
        }

        // ---------------------

        fn evaluate_univariate_term_into_values<F: PrimeField>(
            uni: &UnivariateTerm<F>,
            substituted_witness: &HashMap::<StepDifference<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut HashMap::<WitnessEvaluationData<F>, Polynomial<F, Values>>,
            power_hint: Option<u64>, // make space for later multiplications
            worker: &Worker
        ) -> Result<Polynomial<F, Values>, SynthesisError>
        {
            let scale_to_power = if power_hint.is_some() {
                power_hint.expect("is some") as u64
            } else {
                uni.power
            };

            let base = substituted_witness.get(&uni.steps_difference).expect("should exist").clone();
            let base_len = base.size() as u64;

            let evaluation_data = WitnessEvaluationData {
                mask: uni.steps_difference,
                power: uni.power,
                total_lde_length: scale_to_power*base_len
            };

            if let Some(e) = evaluated_univariate_terms.get(&evaluation_data) {
                return Ok(e.clone());
            }

            let factor = scale_to_power.next_power_of_two() as usize;
            let mut base = base.lde(&worker, factor)?;
            base.pow(&worker, uni.power);

            let one = F::one();
            if uni.coeff != one {
                let mut minus_one = one;
                minus_one.negate();
                if uni.coeff == minus_one {
                    base.negate(&worker);
                } else {
                    base.scale(&worker, uni.coeff);
                }
            }

            evaluated_univariate_terms.insert(evaluation_data, base.clone());

            Ok(base)
        }

        // ---------------------

        // such calls most likely will have start at 0 and num_steps = domain_size - 1
        fn inverse_divisor_for_dense_constraint_in_coset<F: PrimeField> (
            column_domain: &Domain<F>,
            term_evaluation_domain: &Domain<F>,
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
            worker.scope(inverse_divisors.size(), |scope, chunk| {
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

            inverse_divisors.batch_inversion(&worker)?;

            // now do the evaluation

            worker.scope(inverse_divisors.size(), |scope, chunk| {
                for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
                    let roots_iter_outer = roots_iter.clone();
                    scope.spawn(move |_| {
                        let mut x = evaluation_domain_generator.pow([(i*chunk) as u64]);
                        x.mul_assign(&multiplicative_generator);
                        for v in inv_divis.iter_mut() {
                            let mut d = *v;
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

            Ok((inverse_divisors, divisor_degree))
        }

        // ---------------------

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

        // ---------------------

        let worker = Worker::new();
        let num_registers_sup = self.num_registers.next_power_of_two();
        let num_steps_sup = self.num_steps.next_power_of_two();
        let g_size = num_registers_sup * num_steps_sup * self.max_constraint_power;

        let mut g_poly = Polynomial::<F, Coefficients>::new_for_size(g_size)?;
        let subterm_coefficients = Polynomial::<F, Coefficients>::new_for_size(g_size)?;

        let mut current_coeff = F::one();

        let subterm_domain = Domain::new_for_size(subterm_coefficients.size() as u64)?;

        let mut evaluated_terms_map: HashMap::<WitnessEvaluationData<F>, Polynomial<F, Values>> = HashMap::new();

        // one may optimize and save on muptiplications for the most expected case when constraints 
        // all have the same density

        let mut constraints_batched_by_density: HashMap::< ConstraintDensity, Vec<Constraint<F>> > = HashMap::new();

        for constraint in self.constraints.iter() {
            if let Some(batch) = constraints_batched_by_density.get_mut(&constraint.density) {
                batch.push(constraint.clone());
            } else {
                constraints_batched_by_density.insert(constraint.density, vec![constraint.clone()]);
            }
        }

        for (density, constraints) in constraints_batched_by_density.into_iter() {
            let mut per_density_coefficients = subterm_coefficients.clone();

            // TODO: Refactor constraints definitions
            let c0 = constraints[0].clone();
            let start_at = match c0.density {
                ConstraintDensity::Dense => {
                    c0.start_at as u64
                },
                _ => {
                    unimplemented!()
                }
            };
            
            let inverse_divisors = match self.inversed_divisors_in_cosets.get(&density) {
                Some(div) => {
                    div.clone()
                },
                _ => {
                    let (inverse_divisors, _divisor_degree) = match density {
                        ConstraintDensity::Dense => {
                            let result = inverse_divisor_for_dense_constraint_in_coset(
                                &self.column_domain,
                                &subterm_domain, 
                                start_at,
                                self.num_steps as u64,
                                &worker
                            )?;

                            result
                        },
                        _ => {
                            unimplemented!();
                        }
                    };
                    self.inversed_divisors_in_cosets.insert(density, inverse_divisors.clone());

                    inverse_divisors
                }
            };

            for constraint in constraints.into_iter() {
                current_coeff.mul_assign(&alpha);

                for term in constraint.terms.iter() {
                    let evaluated_term = evaluate_constraint_term_into_coefficients(
                        &term, 
                        &self.mask_applied_polynomials,
                        &mut evaluated_terms_map,
                        None,
                        &worker
                    )?;

                    per_density_coefficients.add_assign_scaled(&worker, &evaluated_term, &current_coeff);
                }

                per_density_coefficients.as_mut()[0].add_assign(&constraint.constant_term);
            }

            
            // these values are correct and are evaluations of some polynomial at points (gen, gen * omega, gen * omega*2)
            let mut per_density_in_coset = per_density_coefficients.coset_fft(&worker);
            per_density_in_coset.mul_assign(&worker, &inverse_divisors);

            let per_density_coefficients = per_density_in_coset.icoset_fft(&worker);

            g_poly.add_assign(&worker, &per_density_coefficients);
        }

        // for constraint in self.constraints.iter() {
        //     current_coeff.mul_assign(&alpha);
            
        //     let mut subterm_coefficients = subterm_coefficients.clone();

        //     let mut inverse_divisors = match self.inversed_divisors_in_cosets.get(&constraint.density) {
        //         Some(div) => {
        //             div.clone()
        //         },
        //         _ => {
        //             let (inverse_divisors, _divisor_degree) = match constraint.density {
        //                 ConstraintDensity::Dense => {
        //                     let start_at = constraint.start_at as u64;

        //                     let result = inverse_divisor_for_dense_constraint_in_coset(
        //                         &self.column_domain,
        //                         &subterm_domain, 
        //                         start_at,
        //                         self.num_steps as u64,
        //                         &worker
        //                     )?;

        //                     result
        //                 },
        //                 _ => {
        //                     unimplemented!();
        //                 }
        //             };
        //             self.inversed_divisors_in_cosets.insert(constraint.density, inverse_divisors.clone());

        //             inverse_divisors
        //         }
        //     };

        //     // first we need to calculate denominator at the value

        //     inverse_divisors.scale(&worker, current_coeff);

        //     for term in constraint.terms.iter() {
        //         let evaluated_term = evaluate_constraint_term_into_coefficients(
        //             &term, 
        //             &self.mask_applied_polynomials,
        //             &mut evaluated_terms_map,
        //             None,
        //             &worker
        //         )?;

        //         subterm_coefficients.add_assign(&worker, &evaluated_term);
        //     }

        //     subterm_coefficients.as_mut()[0].add_assign(&constraint.constant_term);

        //     // these values are correct and are evaluations of some polynomial at points (gen, gen * omega, gen * omega*2)
        //     let mut subterm_values_in_coset = subterm_coefficients.coset_fft(&worker);

        //     subterm_values_in_coset.mul_assign(&worker, &inverse_divisors);

        //     let subterm_coefficients = subterm_values_in_coset.icoset_fft(&worker);

        //     g_poly.add_assign(&worker, &subterm_coefficients);
        // }

        for b_constraint in self.boundary_constraints.iter() {
            current_coeff.mul_assign(&alpha);

            // x - a
            let column_generator = self.column_domain.generator;
            let trace_generator = self.full_trace_domain.generator;

            let mut q_poly = Polynomial::<F, Coefficients>::new_for_size(2)?;
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

            let evaluated_term = evaluate_boundary_constraint(
                &b_constraint, 
                &self.mask_applied_polynomials
            )?;

            subterm_coefficients.add_assign(&worker, &evaluated_term);

            let mut inverse_q_poly_coset_values = q_poly.coset_evaluate_at_domain_for_degree_one(
                &worker, 
                subterm_coefficients.size() as u64
            )?;

            inverse_q_poly_coset_values.batch_inversion(&worker)?;
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