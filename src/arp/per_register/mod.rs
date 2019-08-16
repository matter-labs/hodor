use ff::PrimeField;

use super::*;
use crate::domains::Domain;
use crate::SynthesisError;
use crate::polynomials::*;
use crate::fft::multicore::Worker;
use crate::fft::*;
use super::density_query::*;

use super::mappings::*;

fn make_witness_polymonials<F: PrimeField>(
    mut witness: Vec<Vec<F>>, 
    worker: &Worker
) -> Result<Vec<Polynomial<F, Coefficients>>, SynthesisError> {
    let num_registers = witness.len();
    let num_rows = (&witness[0]).len();
    for w in witness.iter() {
        assert!(num_rows == w.len());
    }
    let num_rows_sup = num_rows.next_power_of_two();

    let domain = Domain::<F>::new_for_size(num_rows_sup as u64)?;
    let log_n = domain.power_of_two;
    let omega = domain.generator;
    let omega_inv = omega.inverse().expect("is some");
    let minv = F::from_str(&domain.size.to_string()).expect("exists").inverse().expect("exists");

    let mut result = vec![];

    let num_cpus = worker.cpus;
    let num_cpus_hint = if num_cpus <= num_registers {
        Some(1)
    } else {
        let mut threads_per_register = num_registers / num_cpus;
        if num_registers % num_cpus != 0 {
            threads_per_register += 1;
        }
        Some(threads_per_register)
    };

    worker.scope(0, |scope, _| {
        for w in witness.iter_mut() {
            scope.spawn(move |_| {
                best_fft(&mut w[..], &worker, &omega_inv, log_n as u32, num_cpus_hint);
            });
        }
    });

    for w in witness.into_iter() {
        let mut w = w;
        worker.scope(w.len(), |scope, chunk| {
            for w in w.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for w in w.iter_mut() {
                        w.mul_assign(&minv);
                    }
                });
            }
        });

        let p = Polynomial::from_coeffs(w)?;
        result.push(p);
    }

    Ok(result)
}

impl<F: PrimeField> ARP<F> for ARPInstance<F, PerRegisterARP> {
    fn calculate_witness_polys(
        &self,
        witness: Vec<Vec<F>>,
        worker: &Worker
    ) -> Result<Vec<Polynomial<F, Coefficients>>, SynthesisError> {
        let num_registers = witness.len();
        let num_rows = (&witness[0]).len();
        for w in witness.iter() {
            assert!(num_rows == w.len());
        }

        if num_registers != self.properties.num_registers {
            return Err(SynthesisError::Error);
        }

        if num_rows != self.properties.num_rows {
            return Err(SynthesisError::Error);
        }

        make_witness_polymonials(witness, worker)
    }

    fn from_instance(
        instance: InstanceProperties<F>, 
        _worker: &Worker
    ) -> Result<Self, SynthesisError> {
        let mut t = Self {
            properties: instance,
            _marker: std::marker::PhantomData
        };

        t.route()?;

        Ok(t)
    }

    fn is_satisfied(
        proprties: &InstanceProperties<F>,
        witness: &Vec<Vec<F>>,
        worker: &Worker
    ) -> Result<(), SynthesisError> {
        Self::verify_witness(&proprties, witness, worker)
    }
}

impl<F: PrimeField> ARPInstance<F, PerRegisterARP> {
    /// - make interpolating polynomial f
    /// - add masking coefficients for constraints
    /// - keep boundary constraints as it is
    pub fn route(&mut self) -> Result<(), SynthesisError> {
        let num_rows = self.properties.num_rows as u64;
        let num_rows_sup = num_rows.next_power_of_two();
        let column_domain = Domain::<F>::new_for_size(num_rows_sup)?;

        for mut c in self.properties.constraints.iter_mut() {
            remap_constraint(
                &mut c, 
                &column_domain
            );
        }

        Ok(())
    }

    fn verify_witness(
        properties: &InstanceProperties<F>,
        witness: &Vec<Vec<F>>,
        _worker: &Worker
    ) -> Result<(), SynthesisError> {
        fn evaluate_constraint_on_witness<F: PrimeField>(
            constraint: &Constraint<F>,
            witness: &Vec<Vec<F>>,
            base_row: usize,
        ) -> Result<F, SynthesisError> {
            let mut value = constraint.constant_term;
            for term in constraint.terms.iter() {
                let v = evaluate_term_on_witness(&term, witness, base_row)?;
                value.add_assign(&v);
            }
            
            Ok(value)
        }

        fn evaluate_term_on_witness<F: PrimeField>(
            term: &ConstraintTerm<F>,
            witness: &Vec<Vec<F>>,
            base_row: usize
        ) -> Result<F, SynthesisError> {
            match term {
                ConstraintTerm::Univariate(ref t) => {
                    evaluate_univariate_term_on_witness(
                        t, 
                        &witness,
                        base_row
                    )
                },
                ConstraintTerm::Polyvariate(ref poly_term) => {
                    let mut result = F::one();
                    for t in poly_term.terms.iter() {
                        let v = evaluate_univariate_term_on_witness(
                            &t, 
                            &witness,
                            base_row
                        )?;
                        result.mul_assign(&v);
                    }

                    result.mul_assign(&poly_term.coeff);

                    Ok(result)
                }
            }
        }

        fn evaluate_univariate_term_on_witness<F: PrimeField>(
            univariate_term: &UnivariateTerm<F>,
            witness: &Vec<Vec<F>>,
            base_row: usize
        ) -> Result<F, SynthesisError> {
            let step_delta = match univariate_term.steps_difference {
                StepDifference::Steps(steps) => steps,
                _ => unreachable!()
            };

            let reg_num = match univariate_term.register {
                Register::Register(reg_num) => reg_num,
                _ => unreachable!()
            };

            let num_rows = witness[reg_num].len();
            let row = base_row + step_delta;
            if row >= num_rows {
                return Err(SynthesisError::InvalidValue(
                    format!("access to the value out of the trace at row {}", row)
                    )    
                );
            }

            let value = witness[reg_num][row];
            let mut value = value.pow([univariate_term.power]);
            value.mul_assign(&univariate_term.coeff);

            Ok(value)
        }

        let num_rows = witness[0].len();

        // these constraints are not remapped, so step differences are 
        for c in properties.constraints.iter() {
            let density: Box<dyn DensityQuery> = match c.density {
                ConstraintDensity::Dense(ref dense) => {
                    let dense_query = DenseConstraintQuery::new(&dense, num_rows);

                    Box::from(dense_query)
                },
                _ => {
                    unimplemented!();
                }
            };

            for row in density {
                let value = evaluate_constraint_on_witness(
                    c, 
                    witness, 
                    row
                )?;

                if !value.is_zero() {
                    return Err(SynthesisError::Unsatisfied(
                            format!("constraint:\n{}\nis unsatisfied at the row {}", c, row)
                        )
                    );
                }
            }
        } 

        // these constraints are not remapped, so step differences are 
        for bc in properties.boundary_constraints.iter() {
            let reg_num = match bc.register {
                Register::Register(reg_num) => reg_num,
                _ => unreachable!()
            };

            let at_row = bc.at_row;

            if let Some(expected) = bc.value {
                let from_witness = witness[reg_num][at_row];
                if expected != from_witness {
                    return Err(SynthesisError::Error);
                }
            }
        } 

        Ok(())
    }
}

#[test]
fn test_fib_conversion_into_per_register_arp() {
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::arp::IntoARP;
    use crate::fft::multicore::Worker;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();
    let witness = witness.expect("some witness");
    println!("Witness = {:?}", witness);
    let worker = Worker::new();

    for c in props.constraints.iter() {
        // println!("{:?}", c);
        println!("{}", c);
    }

    let is_satisfied = ARPInstance::<Fr, PerRegisterARP>::is_satisfied(&props, &witness, &worker);
    assert!(is_satisfied.is_ok());
    let arp = ARPInstance::<Fr, PerRegisterARP>::from_instance(props, &worker).expect("must work");

    let witness_polys = arp.calculate_witness_polys(witness.clone(), &worker).expect("must work");
    println!("Witness polys = {:?}", witness_polys);
    for (i, w) in witness_polys.into_iter().enumerate() {
        let vals = w.fft(&worker);
        assert!(witness[i] == vals.into_coeffs());
    }
}

