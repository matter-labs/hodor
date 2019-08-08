use ff::PrimeField;

use super::*;
use crate::domains::Domain;
use crate::SynthesisError;
use crate::polynomials::*;
use crate::fft::multicore::Worker;
use crate::fft::*;

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
        let threads_per_register = (num_registers - 1) / num_cpus + 1;
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
}

impl<F: PrimeField> ARPInstance<F, PerRegisterARP> {
    /// - make interpolating polynomial f
    /// - add masking coefficients for constraints
    /// - keep boundary constraints as it is
    pub fn route(&mut self) -> Result<(), SynthesisError> {
        let num_rows = self.properties.num_rows as u64;
        let num_rows_sup = num_rows.next_power_of_two();
        let column_domain = Domain::<F>::new_for_size(num_rows_sup)?;

        fn remap_univariate_term<F: PrimeField>(
            term: &mut UnivariateTerm<F>,
            column_domain: &Domain<F>,
        ) {
            let step_delta = match term.steps_difference {
                StepDifference::Steps(num) => {
                    num as u64
                },
                _ => {
                    unreachable!("Step differences are not masks yet");
                }
            };

            let mask = column_domain.generator.pow([step_delta]);

            term.steps_difference = StepDifference::Mask(mask);
        }

        fn remap_term<F: PrimeField>(
            term: &mut ConstraintTerm<F>,
            column_domain: &Domain<F>,
        ) {
            match term {
                ConstraintTerm::Univariate(ref mut t) => {
                    remap_univariate_term(
                        t, 
                        &column_domain
                    );
                },
                ConstraintTerm::Polyvariate(ref mut poly_term) => {
                    for mut t in poly_term.terms.iter_mut() {
                        remap_univariate_term(
                            &mut t, 
                            &column_domain
                        );
                    }
                }
            }
        }

        fn remap_constraint<F: PrimeField>(
            constraint: &mut Constraint<F>,
            column_domain: &Domain<F>,
        ) {
            for mut t in constraint.terms.iter_mut() {
                remap_term(
                    &mut t, 
                    &column_domain
                );
            }
        }

        for mut c in self.properties.constraints.iter_mut() {
            remap_constraint(
                &mut c, 
                &column_domain
            );
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
    let arp = ARPInstance::<Fr, PerRegisterARP>::from_instance(props, &worker).expect("must work");
    for c in arp.properties.constraints.iter() {
        println!("{:?}", c);
    }

    let witness_polys = arp.calculate_witness_polys(witness.clone(), &worker).expect("must work");
    println!("Witness polys = {:?}", witness_polys);
    for (i, w) in witness_polys.into_iter().enumerate() {
        let vals = w.fft(&worker);
        assert!(witness[i] == vals.into_coeffs());
    }
}

