use ff::PrimeField;

use crate::air::*;
use crate::arp::*;
use crate::domains::*;
use crate::fft::cooley_tukey_ntt::*;
use crate::fft::multicore::Worker;
use crate::polynomials::*;
use crate::precomputations::*;
use crate::transcript::Transcript;
use crate::SynthesisError;

use indexmap::IndexMap;
use indexmap::IndexSet;

// use std::collections::{IndexMap, IndexSet};

use super::*;

pub mod deep;

pub struct ALIInstance<F: PrimeField, T: ARPType> {
    pub properties: InstanceProperties<F>,
    pub max_constraint_power: u64,
    pub column_domain: Domain<F>,
    pub constraints_domain: Domain<F>,
    pub(crate) all_masks: IndexSet<MaskProperties<F>>,
    pub all_boundary_constrained_registers: IndexSet<Register>,
    pub constraint_divisors: IndexMap<Box<dyn ConstraintDensity<F>>, Polynomial<F, Values>>,
    pub boundary_constraint_divisors: IndexMap<u64, Polynomial<F, Values>>,
    pub constraints_batched_by_density: IndexMap<Box<dyn ConstraintDensity<F>>, Vec<Constraint<F>>>,
    pub precomputations: PrecomputedOmegas<F>,
    _marker: std::marker::PhantomData<T>,
}

impl<F: PrimeField> ALIInstance<F, PerRegisterARP> {
    pub fn from_arp(
        arp: ARPInstance<F, PerRegisterARP>,
        worker: &Worker,
    ) -> Result<ALIInstance<F, PerRegisterARP>, SynthesisError> {
        let mut max_constraint_power: u64 = 0;
        for c in arp.properties.constraints.iter() {
            if c.degree > max_constraint_power {
                max_constraint_power = c.degree;
            }
        }

        let num_rows = arp.properties.num_rows;

        let column_domain = Domain::<F>::new_for_size(arp.properties.num_rows as u64)?;
        let constraints_domain =
            Domain::<F>::new_for_size(column_domain.size * max_constraint_power)?;
        let precomputed_omegas = PrecomputedOmegas::new_for_domain(&constraints_domain, &worker);
        let mut all_masks: IndexSet<MaskProperties<F>> = IndexSet::new();

        for c in arp.properties.constraints.iter() {
            get_masks_from_constraint(&mut all_masks, c);
        }

        // such calls most likely will have start at 0 and num_steps = domain_size - 1
        fn inverse_divisor_for_dense_constraint_in_coset<F: PrimeField>(
            column_domain: &Domain<F>,
            evaluation_domain: &Domain<F>,
            dense_constraint: DenseConstraint,
            // precomputations: PrecomputedOmegas<F>,
            num_rows: u64,
            worker: &Worker,
        ) -> Result<(Polynomial<F, Values>, usize), SynthesisError> {
            let start_at = dense_constraint.start_at;
            let span = dense_constraint.span as u64;
            let mut divisor_degree = column_domain.size as usize;
            let divisor_domain_size = column_domain.size;
            divisor_degree -= start_at as usize;
            divisor_degree -= (divisor_domain_size - num_rows) as usize;
            divisor_degree -= span as usize;

            let roots = {
                let roots_generator = column_domain.generator;

                let mut roots = vec![];
                let mut root = F::one();
                for _ in 0..start_at {
                    roots.push(root);
                    root.mul_assign(&roots_generator);
                }

                let last_step = num_rows - span;
                let mut root = roots_generator.pow([last_step]);
                for _ in last_step..divisor_domain_size {
                    roots.push(root);
                    root.mul_assign(&roots_generator);
                }

                roots
            };

            let roots_iter = roots.iter();

            let evaluation_domain_generator = evaluation_domain.generator;
            let multiplicative_generator = F::multiplicative_generator();

            let evaluation_domain_size = evaluation_domain.size as usize;
            // assert!(evaluation_domain_size == precomputations.coset.len());

            // strategy:
            // - pick evaluation domain (will depend on constraint power) and column domain
            // - evaluate polynomial {X}^T - 1 where T is a size of the column domain over the coset of evaluation domain
            // - now we need to multiply those values by ({X} - root) for roots of this constraint
            // - first inverse {X}^T - 1 values
            // - for each of these values evaluate ({X} - root) and muptiply by this value

            // TODO: check if we need precomputation at all

            // these are values at the coset
            let mut inverse_divisors =
                Polynomial::<F, Values>::new_for_size(evaluation_domain_size)?;

            // prepare for batch inversion
            worker.scope(inverse_divisors.size(), |scope, chunk| {
                for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
                    scope.spawn(move |_| {
                        let mut x = evaluation_domain_generator.pow([(i * chunk) as u64]);
                        x.mul_assign(&multiplicative_generator);
                        for v in inv_divis.iter_mut() {
                            *v = x.pow([divisor_domain_size]);
                            v.sub_assign(&F::one());

                            x.mul_assign(&evaluation_domain_generator);
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
                        let mut x = evaluation_domain_generator.pow([(i * chunk) as u64]);
                        x.mul_assign(&multiplicative_generator);
                        for v in inv_divis.iter_mut() {
                            let mut d = *v;
                            for root in roots_iter_outer.clone() {
                                // (X - root)
                                let mut tmp = x;
                                tmp.sub_assign(&root);
                                d.mul_assign(&tmp);
                            }
                            // 1 / ( (X^T-1) / (X - 1)(X - omega)(...) ) =  (X - 1)(X - omega)(...) / (X^T-1)
                            *v = d;

                            x.mul_assign(&evaluation_domain_generator);
                        }
                    });
                }
            });

            Ok((inverse_divisors, divisor_degree))
        }

        let mut constraints_batched_by_density: IndexMap<
            Box<dyn ConstraintDensity<F>>,
            Vec<Constraint<F>>,
        > = IndexMap::new();

        for constraint in arp.properties.constraints.iter() {
            if let Some(batch) = constraints_batched_by_density.get_mut(&constraint.density) {
                batch.push(constraint.clone());
            } else {
                constraints_batched_by_density
                    .insert(constraint.density.clone(), vec![constraint.clone()]);
            }
        }

        let mut inverse_divisors_per_constraint_density: IndexMap<
            Box<dyn ConstraintDensity<F>>,
            Polynomial<F, Values>,
        > = IndexMap::new();

        for (density, _) in constraints_batched_by_density.iter() {
            let (inverse_divisors, _) = density.inverse_divisor_in_coset(
                &column_domain,
                &constraints_domain,
                &None,
                &None,
                num_rows,
                &worker,
            )?;

            inverse_divisors_per_constraint_density.insert(density.clone(), inverse_divisors);
        }

        let mut all_boundary_constrained_registers: IndexSet<Register> = IndexSet::new();

        for b_c in arp.properties.boundary_constraints.iter() {
            if all_boundary_constrained_registers
                .get(&b_c.register)
                .is_none()
            {
                all_boundary_constrained_registers.insert(b_c.register);
            }
        }

        let mut boundary_constraint_divisors: IndexMap<u64, Polynomial<F, Values>> =
            IndexMap::new();

        let mut all_boundary_constraint_steps: IndexSet<usize> = IndexSet::new();

        for b_c in arp.properties.boundary_constraints.iter() {
            get_mask_from_boundary_constraint(&mut all_masks, b_c);

            if all_boundary_constraint_steps.get(&b_c.at_row).is_none() {
                all_boundary_constraint_steps.insert(b_c.at_row);
            }
        }

        for row in all_boundary_constraint_steps.into_iter() {
            // precompute divisors
            let mut q_poly = Polynomial::<F, Coefficients>::new_for_size(2)?;
            q_poly.as_mut()[1] = F::one();
            let root = column_domain.generator.pow([row as u64]);
            q_poly.as_mut()[0].sub_assign(&root);
            let mut inverse_q_poly_coset_values = q_poly
                .coset_evaluate_at_domain_for_degree_one(&worker, constraints_domain.size as u64)?;

            inverse_q_poly_coset_values.batch_inversion(&worker)?;
            boundary_constraint_divisors.insert(row as u64, inverse_q_poly_coset_values);
        }

        let result = ALIInstance::<F, PerRegisterARP> {
            properties: arp.properties,
            max_constraint_power: max_constraint_power,
            column_domain: column_domain,
            constraints_domain: constraints_domain,
            all_masks: all_masks,
            all_boundary_constrained_registers: all_boundary_constrained_registers,
            constraint_divisors: inverse_divisors_per_constraint_density,
            boundary_constraint_divisors: boundary_constraint_divisors,
            constraints_batched_by_density: constraints_batched_by_density,
            precomputations: precomputed_omegas,
            _marker: std::marker::PhantomData,
        };

        Ok(result)
    }

    pub fn calculate_g<T: Transcript<F>>(
        &self,
        transcript: &mut T,
        witness: Vec<Polynomial<F, Coefficients>>,
        worker: &Worker,
    ) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
        assert!(witness.len() == self.properties.num_registers);

        // all zeroes
        let mut g_values =
            Polynomial::<F, Values>::new_for_size(self.constraints_domain.size as usize)?;

        fn evaluate_for_mask<F: PrimeField>(
            mut f: Polynomial<F, Coefficients>,
            mask: StepDifference<F>,
            worker: &Worker,
        ) -> Polynomial<F, Coefficients> {
            match mask {
                StepDifference::Mask(m) => {
                    f.distribute_powers(&worker, m);
                }
                _ => {
                    unreachable!();
                }
            }

            f
        }

        // first mask all the registers

        let mut mask_applied_polynomials =
            IndexMap::<MaskProperties<F>, Polynomial<F, Coefficients>, _>::new();

        for m in self.all_masks.iter() {
            let register_number = match m.register {
                Register::Register(n) => n,
                _ => {
                    unreachable!();
                }
            };
            let poly = (&witness[register_number]).clone();
            let masked = evaluate_for_mask(poly, m.steps_difference, &worker);
            mask_applied_polynomials.insert(m.clone(), masked);
        }

        fn calculate_adjustment_polynomial_in_coset<F: PrimeField>(
            adjustment: u64,
            alpha: F,
            beta: F,
            domain: &Domain<F>,
            precomputations: &PrecomputedOmegas<F>,
            worker: &Worker,
        ) -> Polynomial<F, Values> {
            assert!(adjustment >= 1);
            assert!(precomputations.coset.len() as u64 == domain.size);
            let mut poly = Polynomial::from_values(precomputations.coset.clone()).expect("is ok");
            poly.pow(&worker, adjustment);
            poly.scale(&worker, alpha);
            poly.add_constant(&worker, &beta);

            poly
        }

        // returns constraint evaluated in the coset
        fn evaluate_constraint_term_into_values<F: PrimeField>(
            term: &ConstraintTerm<F>,
            substituted_witness: &IndexMap<MaskProperties<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut IndexMap<
                WitnessEvaluationData<F>,
                Polynomial<F, Values>,
            >,
            power_hint: u64,
            worker: &Worker,
        ) -> Result<Polynomial<F, Values>, SynthesisError> {
            assert!(power_hint.is_power_of_two());
            let result = match term {
                ConstraintTerm::Univariate(uni) => {
                    let t = evaluate_univariate_term_into_values(
                        uni,
                        substituted_witness,
                        evaluated_univariate_terms,
                        power_hint,
                        worker,
                    )?;

                    t
                }
                ConstraintTerm::Polyvariate(poly) => {
                    let mut values_result: Option<Polynomial<F, Values>> = None;
                    // evaluate subcomponents in a value form and multiply
                    for uni in poly.terms.iter() {
                        let t = evaluate_univariate_term_into_values(
                            uni,
                            substituted_witness,
                            evaluated_univariate_terms,
                            power_hint,
                            &worker,
                        )?;
                        if let Some(res) = values_result.as_mut() {
                            res.mul_assign(&worker, &t);
                        } else {
                            values_result = Some(t);
                        }
                    }

                    let mut as_values = values_result.expect("is some");
                    as_values.scale(&worker, poly.coeff);

                    as_values
                }
            };

            Ok(result)
        }

        // ---------------------

        // returns univariate term evaluated at coset
        fn evaluate_univariate_term_into_values<F: PrimeField>(
            uni: &UnivariateTerm<F>,
            substituted_witness: &IndexMap<MaskProperties<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut IndexMap<
                WitnessEvaluationData<F>,
                Polynomial<F, Values>,
            >,
            power_hint: u64,
            worker: &Worker,
        ) -> Result<Polynomial<F, Values>, SynthesisError> {
            assert!(power_hint.is_power_of_two());
            let mask_props = MaskProperties::<F> {
                register: uni.register,
                steps_difference: uni.steps_difference,
            };
            let base = substituted_witness
                .get(&mask_props)
                .expect("should exist")
                .clone();
            let base_len = base.size() as u64;

            let evaluation_data = WitnessEvaluationData::<F> {
                mask: mask_props,
                power: uni.power,
                total_lde_length: power_hint * base_len,
            };

            if let Some(e) = evaluated_univariate_terms.get(&evaluation_data) {
                let mut base = e.clone();
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
                return Ok(base);
            }

            let factor = power_hint as usize;
            assert!(factor.is_power_of_two());
            let mut base = base.coset_lde(&worker, factor)?;
            // let mut base = base.lde(&worker, factor)?;
            base.pow(&worker, uni.power);

            evaluated_univariate_terms.insert(evaluation_data, base.clone());

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

            Ok(base)
        }

        let mut evaluated_terms_map: IndexMap<WitnessEvaluationData<F>, Polynomial<F, Values>> =
            IndexMap::new();

        // now evaluate TF constraints
        for (density, batch) in self.constraints_batched_by_density.iter() {
            let mut batch_values = g_values.clone();
            for c in batch.iter() {
                let constraint_power = c.degree;
                assert!(self.max_constraint_power >= constraint_power);
                let adjustment = self.max_constraint_power - constraint_power;
                let alpha = transcript.get_challenge();
                let beta = transcript.get_challenge();

                let adj_poly = if adjustment == 0 {
                    None
                } else {
                    let adj_poly = calculate_adjustment_polynomial_in_coset(
                        adjustment,
                        alpha,
                        beta,
                        &self.constraints_domain,
                        &self.precomputations,
                        &worker,
                    );

                    Some(adj_poly)
                };

                let mut constraint_values = g_values.clone();
                for t in c.terms.iter() {
                    let subval = evaluate_constraint_term_into_values(
                        t,
                        &mask_applied_polynomials,
                        &mut evaluated_terms_map,
                        self.max_constraint_power,
                        &worker,
                    )?;
                    constraint_values.add_assign(&worker, &subval);
                }

                constraint_values.add_constant(&worker, &c.constant_term);
                if let Some(adj) = adj_poly {
                    constraint_values.mul_assign(&worker, &adj);
                } else {
                    // just apply alpha
                    constraint_values.scale(&worker, alpha);
                }

                batch_values.add_assign(&worker, &constraint_values);
            }

            let divisors = self.constraint_divisors.get(density).expect("is some");

            batch_values.mul_assign(&worker, divisors);

            g_values.add_assign(&worker, &batch_values);
        }

        let boundary_lde_factor = self.max_constraint_power;

        // now evaluate normal constraints
        for b_c in self.properties.boundary_constraints.iter() {
            let alpha = transcript.get_challenge();
            let beta = transcript.get_challenge();
            let adjustment = self.max_constraint_power - 1;

            let adj_poly = if adjustment == 0 {
                None
            } else {
                let adj_poly = calculate_adjustment_polynomial_in_coset(
                    adjustment,
                    alpha,
                    beta,
                    &self.constraints_domain,
                    &self.precomputations,
                    &worker,
                );

                Some(adj_poly)
            };

            let reg_num = match b_c.register {
                Register::Register(reg_number) => reg_number,
                _ => {
                    unreachable!();
                }
            };
            let mut witness_poly = (&witness[reg_num]).clone();
            witness_poly.as_mut()[0].sub_assign(&b_c.value.expect("is some"));
            let mut constraint_values =
                witness_poly.coset_lde(&worker, boundary_lde_factor as usize)?;
            if let Some(adj) = adj_poly {
                constraint_values.mul_assign(&worker, &adj);
            } else {
                // just apply alpha
                constraint_values.scale(&worker, alpha);
            }

            let divisors = self
                .boundary_constraint_divisors
                .get(&(b_c.at_row as u64))
                .expect("is some");
            constraint_values.mul_assign(&worker, divisors);

            g_values.add_assign(&worker, &constraint_values);
        }

        let g_poly = g_values.icoset_fft(&worker);

        Ok(g_poly)
    }
    pub fn calculate_g_with_ntt<T: Transcript<F>>(
        &self,
        transcript: &mut T,
        witness: Vec<Polynomial<F, Coefficients>>,
        worker: &Worker,
        bitreversed_omegas: &BitReversedOmegas<F>,
        coset_factor: &F,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        assert!(witness.len() == self.properties.num_registers);

        // all zeroes
        let mut g_values =
            Polynomial::<F, Values>::new_for_size(self.constraints_domain.size as usize)?;

        fn evaluate_for_mask<F: PrimeField>(
            mut f: Polynomial<F, Coefficients>,
            mask: StepDifference<F>,
            worker: &Worker,
        ) -> Polynomial<F, Coefficients> {
            match mask {
                StepDifference::Mask(m) => {
                    f.distribute_powers(&worker, m);
                }
                _ => {
                    unreachable!();
                }
            }

            f
        }

        // first mask all the registers

        let mut mask_applied_polynomials =
            IndexMap::<MaskProperties<F>, Polynomial<F, Coefficients>, _>::new();

        for m in self.all_masks.iter() {
            let register_number = match m.register {
                Register::Register(n) => n,
                _ => {
                    unreachable!();
                }
            };
            let poly = (&witness[register_number]).clone();
            let masked = evaluate_for_mask(poly, m.steps_difference, &worker);
            mask_applied_polynomials.insert(m.clone(), masked);
        }

        fn calculate_adjustment_polynomial_in_coset<F: PrimeField>(
            adjustment: u64,
            alpha: F,
            beta: F,
            domain: &Domain<F>,
            precomputations: &PrecomputedOmegas<F>,
            worker: &Worker,
        ) -> Polynomial<F, Values> {
            assert!(adjustment >= 1);
            assert!(precomputations.coset.len() as u64 == domain.size);
            let mut poly = Polynomial::from_values(precomputations.coset.clone()).expect("is ok");
            poly.pow(&worker, adjustment);
            poly.scale(&worker, alpha);
            poly.add_constant(&worker, &beta);

            poly
        }

        // returns constraint evaluated in the coset
        fn evaluate_constraint_term_into_values<F: PrimeField>(
            term: &ConstraintTerm<F>,
            substituted_witness: &IndexMap<MaskProperties<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut IndexMap<
                WitnessEvaluationData<F>,
                Polynomial<F, Values>,
            >,
            power_hint: u64,
            worker: &Worker,
        ) -> Result<Polynomial<F, Values>, SynthesisError> {
            assert!(power_hint.is_power_of_two());
            let result = match term {
                ConstraintTerm::Univariate(uni) => {
                    let t = evaluate_univariate_term_into_values(
                        uni,
                        substituted_witness,
                        evaluated_univariate_terms,
                        power_hint,
                        worker,
                    )?;

                    t
                }
                ConstraintTerm::Polyvariate(poly) => {
                    let mut values_result: Option<Polynomial<F, Values>> = None;
                    // evaluate subcomponents in a value form and multiply
                    for uni in poly.terms.iter() {
                        let t = evaluate_univariate_term_into_values(
                            uni,
                            substituted_witness,
                            evaluated_univariate_terms,
                            power_hint,
                            &worker,
                        )?;
                        if let Some(res) = values_result.as_mut() {
                            res.mul_assign(&worker, &t);
                        } else {
                            values_result = Some(t);
                        }
                    }

                    let mut as_values = values_result.expect("is some");
                    as_values.scale(&worker, poly.coeff);

                    as_values
                }
            };

            Ok(result)
        }

        // ---------------------

        // returns univariate term evaluated at coset
        fn evaluate_univariate_term_into_values<F: PrimeField>(
            uni: &UnivariateTerm<F>,
            substituted_witness: &IndexMap<MaskProperties<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut IndexMap<
                WitnessEvaluationData<F>,
                Polynomial<F, Values>,
            >,
            power_hint: u64,
            worker: &Worker,
        ) -> Result<Polynomial<F, Values>, SynthesisError> {
            assert!(power_hint.is_power_of_two());
            let mask_props = MaskProperties::<F> {
                register: uni.register,
                steps_difference: uni.steps_difference,
            };
            let base = substituted_witness
                .get(&mask_props)
                .expect("should exist")
                .clone();
            let base_len = base.size() as u64;

            let evaluation_data = WitnessEvaluationData::<F> {
                mask: mask_props,
                power: uni.power,
                total_lde_length: power_hint * base_len,
            };

            if let Some(e) = evaluated_univariate_terms.get(&evaluation_data) {
                let mut base = e.clone();
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
                return Ok(base);
            }

            let factor = power_hint as usize;
            assert!(factor.is_power_of_two());
            let mut base = base.coset_lde(&worker, factor)?;
            // let mut base = base.lde(&worker, factor)?;
            base.pow(&worker, uni.power);

            evaluated_univariate_terms.insert(evaluation_data, base.clone());

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

            Ok(base)
        }

        let mut evaluated_terms_map: IndexMap<WitnessEvaluationData<F>, Polynomial<F, Values>> =
            IndexMap::new();

        // now evaluate TF constraints
        for (density, batch) in self.constraints_batched_by_density.iter() {
            let mut batch_values = g_values.clone();
            for c in batch.iter() {
                let constraint_power = c.degree;
                assert!(self.max_constraint_power >= constraint_power);
                let adjustment = self.max_constraint_power - constraint_power;
                let alpha = transcript.get_challenge();
                let beta = transcript.get_challenge();

                let adj_poly = if adjustment == 0 {
                    None
                } else {
                    let adj_poly = calculate_adjustment_polynomial_in_coset(
                        adjustment,
                        alpha,
                        beta,
                        &self.constraints_domain,
                        &self.precomputations,
                        &worker,
                    );

                    Some(adj_poly)
                };

                let mut constraint_values = g_values.clone();
                for t in c.terms.iter() {
                    let subval = evaluate_constraint_term_into_values(
                        t,
                        &mask_applied_polynomials,
                        &mut evaluated_terms_map,
                        self.max_constraint_power,
                        &worker,
                    )?;
                    constraint_values.add_assign(&worker, &subval);
                }

                constraint_values.add_constant(&worker, &c.constant_term);
                if let Some(adj) = adj_poly {
                    constraint_values.mul_assign(&worker, &adj);
                } else {
                    // just apply alpha
                    constraint_values.scale(&worker, alpha);
                }

                batch_values.add_assign(&worker, &constraint_values);
            }

            let divisors = self.constraint_divisors.get(density).expect("is some");

            batch_values.mul_assign(&worker, divisors);

            g_values.add_assign(&worker, &batch_values);
        }

        let boundary_lde_factor = self.max_constraint_power;

        // now evaluate normal constraints
        for b_c in self.properties.boundary_constraints.iter() {
            let alpha = transcript.get_challenge();
            let beta = transcript.get_challenge();
            let adjustment = self.max_constraint_power - 1;

            let adj_poly = if adjustment == 0 {
                None
            } else {
                let adj_poly = calculate_adjustment_polynomial_in_coset(
                    adjustment,
                    alpha,
                    beta,
                    &self.constraints_domain,
                    &self.precomputations,
                    &worker,
                );

                Some(adj_poly)
            };

            let reg_num = match b_c.register {
                Register::Register(reg_number) => reg_number,
                _ => {
                    unreachable!();
                }
            };
            let mut witness_poly = (&witness[reg_num]).clone();
            witness_poly.as_mut()[0].sub_assign(&b_c.value.expect("is some"));
            let mut constraint_values = witness_poly.bitreversed_lde_using_bitreversed_ntt(
                &worker,
                boundary_lde_factor as usize,
                bitreversed_omegas,
                coset_factor,
            )?;
            if let Some(adj) = adj_poly {
                constraint_values.mul_assign(&worker, &adj);
            } else {
                // just apply alpha
                constraint_values.scale(&worker, alpha);
            }

            let divisors = self
                .boundary_constraint_divisors
                .get(&(b_c.at_row as u64))
                .expect("is some");
            constraint_values.mul_assign(&worker, divisors);

            g_values.add_assign(&worker, &constraint_values);
        }

        // let g_poly = g_values.icoset_fft(&worker);
        

        Ok(g_values)
    }
    pub fn calculate_g_with_square_root_fft<T: Transcript<F>>(
        &self,
        transcript: &mut T,
        witness: Vec<Polynomial<F, Coefficients>>,
        worker: &Worker,
        precomputed_twiddle_factors: &Vec<F>,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        assert!(witness.len() == self.properties.num_registers);

        // all zeroes
        let mut g_values =
            Polynomial::<F, Values>::new_for_size(self.constraints_domain.size as usize)?;

        fn evaluate_for_mask<F: PrimeField>(
            mut f: Polynomial<F, Coefficients>,
            mask: StepDifference<F>,
            worker: &Worker,
        ) -> Polynomial<F, Coefficients> {
            match mask {
                StepDifference::Mask(m) => {
                    f.distribute_powers(&worker, m);
                }
                _ => {
                    unreachable!();
                }
            }

            f
        }

        // first mask all the registers

        let mut mask_applied_polynomials =
            IndexMap::<MaskProperties<F>, Polynomial<F, Coefficients>, _>::new();

        for m in self.all_masks.iter() {
            let register_number = match m.register {
                Register::Register(n) => n,
                _ => {
                    unreachable!();
                }
            };
            let poly = (&witness[register_number]).clone();
            let masked = evaluate_for_mask(poly, m.steps_difference, &worker);
            mask_applied_polynomials.insert(m.clone(), masked);
        }

        fn calculate_adjustment_polynomial_in_coset<F: PrimeField>(
            adjustment: u64,
            alpha: F,
            beta: F,
            domain: &Domain<F>,
            precomputations: &PrecomputedOmegas<F>,
            worker: &Worker,
        ) -> Polynomial<F, Values> {
            assert!(adjustment >= 1);
            assert!(precomputations.coset.len() as u64 == domain.size);
            let mut poly = Polynomial::from_values(precomputations.coset.clone()).expect("is ok");
            poly.pow(&worker, adjustment);
            poly.scale(&worker, alpha);
            poly.add_constant(&worker, &beta);

            poly
        }

        // returns constraint evaluated in the coset
        fn evaluate_constraint_term_into_values<F: PrimeField>(
            term: &ConstraintTerm<F>,
            substituted_witness: &IndexMap<MaskProperties<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut IndexMap<
                WitnessEvaluationData<F>,
                Polynomial<F, Values>,
            >,
            power_hint: u64,
            worker: &Worker,
        ) -> Result<Polynomial<F, Values>, SynthesisError> {
            assert!(power_hint.is_power_of_two());
            let result = match term {
                ConstraintTerm::Univariate(uni) => {
                    let t = evaluate_univariate_term_into_values(
                        uni,
                        substituted_witness,
                        evaluated_univariate_terms,
                        power_hint,
                        worker,
                    )?;

                    t
                }
                ConstraintTerm::Polyvariate(poly) => {
                    let mut values_result: Option<Polynomial<F, Values>> = None;
                    // evaluate subcomponents in a value form and multiply
                    for uni in poly.terms.iter() {
                        let t = evaluate_univariate_term_into_values(
                            uni,
                            substituted_witness,
                            evaluated_univariate_terms,
                            power_hint,
                            &worker,
                        )?;
                        if let Some(res) = values_result.as_mut() {
                            res.mul_assign(&worker, &t);
                        } else {
                            values_result = Some(t);
                        }
                    }

                    let mut as_values = values_result.expect("is some");
                    as_values.scale(&worker, poly.coeff);

                    as_values
                }
            };

            Ok(result)
        }

        // ---------------------

        // returns univariate term evaluated at coset
        fn evaluate_univariate_term_into_values<F: PrimeField>(
            uni: &UnivariateTerm<F>,
            substituted_witness: &IndexMap<MaskProperties<F>, Polynomial<F, Coefficients>>,
            evaluated_univariate_terms: &mut IndexMap<
                WitnessEvaluationData<F>,
                Polynomial<F, Values>,
            >,
            power_hint: u64,
            worker: &Worker,
        ) -> Result<Polynomial<F, Values>, SynthesisError> {
            assert!(power_hint.is_power_of_two());
            let mask_props = MaskProperties::<F> {
                register: uni.register,
                steps_difference: uni.steps_difference,
            };
            let base = substituted_witness
                .get(&mask_props)
                .expect("should exist")
                .clone();
            let base_len = base.size() as u64;

            let evaluation_data = WitnessEvaluationData::<F> {
                mask: mask_props,
                power: uni.power,
                total_lde_length: power_hint * base_len,
            };

            if let Some(e) = evaluated_univariate_terms.get(&evaluation_data) {
                let mut base = e.clone();
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
                return Ok(base);
            }

            let factor = power_hint as usize;
            assert!(factor.is_power_of_two());
            let mut base = base.coset_lde(&worker, factor)?;
            // let mut base = base.lde(&worker, factor)?;
            base.pow(&worker, uni.power);

            evaluated_univariate_terms.insert(evaluation_data, base.clone());

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

            Ok(base)
        }

        let mut evaluated_terms_map: IndexMap<WitnessEvaluationData<F>, Polynomial<F, Values>> =
            IndexMap::new();

        // now evaluate TF constraints
        for (density, batch) in self.constraints_batched_by_density.iter() {
            let mut batch_values = g_values.clone();
            for c in batch.iter() {
                let constraint_power = c.degree;
                assert!(self.max_constraint_power >= constraint_power);
                let adjustment = self.max_constraint_power - constraint_power;
                let alpha = transcript.get_challenge();
                let beta = transcript.get_challenge();

                let adj_poly = if adjustment == 0 {
                    None
                } else {
                    let adj_poly = calculate_adjustment_polynomial_in_coset(
                        adjustment,
                        alpha,
                        beta,
                        &self.constraints_domain,
                        &self.precomputations,
                        &worker,
                    );

                    Some(adj_poly)
                };

                let mut constraint_values = g_values.clone();
                for t in c.terms.iter() {
                    let subval = evaluate_constraint_term_into_values(
                        t,
                        &mask_applied_polynomials,
                        &mut evaluated_terms_map,
                        self.max_constraint_power,
                        &worker,
                    )?;
                    constraint_values.add_assign(&worker, &subval);
                }

                constraint_values.add_constant(&worker, &c.constant_term);
                if let Some(adj) = adj_poly {
                    constraint_values.mul_assign(&worker, &adj);
                } else {
                    // just apply alpha
                    constraint_values.scale(&worker, alpha);
                }

                batch_values.add_assign(&worker, &constraint_values);
            }

            let divisors = self.constraint_divisors.get(density).expect("is some");

            batch_values.mul_assign(&worker, divisors);

            g_values.add_assign(&worker, &batch_values);
        }

        let boundary_lde_factor = self.max_constraint_power;

        // now evaluate normal constraints
        for b_c in self.properties.boundary_constraints.iter() {
            let alpha = transcript.get_challenge();
            let beta = transcript.get_challenge();
            let adjustment = self.max_constraint_power - 1;

            let adj_poly = if adjustment == 0 {
                None
            } else {
                let adj_poly = calculate_adjustment_polynomial_in_coset(
                    adjustment,
                    alpha,
                    beta,
                    &self.constraints_domain,
                    &self.precomputations,
                    &worker,
                );

                Some(adj_poly)
            };

            let reg_num = match b_c.register {
                Register::Register(reg_number) => reg_number,
                _ => {
                    unreachable!();
                }
            };
            let mut witness_poly = (&witness[reg_num]).clone();
            witness_poly.as_mut()[0].sub_assign(&b_c.value.expect("is some"));
            let mut constraint_values = witness_poly.lde_using_square_root_fft(
                &worker,
                boundary_lde_factor as usize,
                precomputed_twiddle_factors,
            )?;
            if let Some(adj) = adj_poly {
                constraint_values.mul_assign(&worker, &adj);
            } else {
                // just apply alpha
                constraint_values.scale(&worker, alpha);
            }

            let divisors = self
                .boundary_constraint_divisors
                .get(&(b_c.at_row as u64))
                .expect("is some");
            constraint_values.mul_assign(&worker, divisors);

            g_values.add_assign(&worker, &constraint_values);
        }        

        Ok(g_values)
    }

}

#[test]
fn test_fib_conversion_into_ali() {
    use crate::air::Fibonacci;
    use crate::air::IntoAIR;
    use crate::air::TestTraceSystem;
    use crate::arp::IntoARP;
    use crate::fft::multicore::Worker;
    use crate::transcript::Transcript;
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
    let witness = witness.expect("some witness");
    // println!("Witness = {:?}", witness);
    let worker = Worker::new();
    let arp = ARPInstance::<Fr, PerRegisterARP>::from_instance(props, &worker).expect("must work");
    let witness_polys = arp
        .calculate_witness_polys(witness, &worker)
        .expect("must work");
    // println!("Witness polys = {:?}", witness_polys);

    let ali = ALIInstance::from_arp(arp, &worker).expect("is some");
    let mut transcript = crate::transcript::Blake2sTranscript::new();
    let alpha = Fr::from_str("123").unwrap();
    transcript.commit_field_element(&alpha);

    let g_poly_interpolant = ali
        .calculate_g(&mut transcript, witness_polys, &worker)
        .expect("is some");

    println!("G coefficients = {:?}", g_poly_interpolant);
    assert!(g_poly_interpolant.as_ref()[3].is_zero());
    let g_values = g_poly_interpolant.fft(&worker);
    println!("G values = {:?}", g_values);
}
