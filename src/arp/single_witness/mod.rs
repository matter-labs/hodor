use crate::air::*;
use ff::PrimeField;

use super::*;
use crate::domains::Domain;
use crate::SynthesisError;
use crate::polynomials::*;
use crate::fft::multicore::Worker;

// ARP works with remapped registers and no longer cares about their meaning
#[derive(Clone, Debug)]
pub struct SingleWitnessARPInstance<F:PrimeField> {
    pub witness_poly: Option<Polynomial<F, Coefficients>>,
    pub num_rows: usize,
    pub num_registers: usize,
    pub constraints: Vec<Constraint<F>>,
    pub boundary_constraints: Vec<BoundaryConstraint<F>>
}

fn make_single_witness_polymonial<F: PrimeField>(
    witness: Vec<Vec<F>>, 
    worker: &Worker
) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    let num_rows = (&witness[0]).len();
    let num_registers = witness.len();

    let num_registers_sup = num_registers.next_power_of_two();
    let num_rows_sup = num_rows.next_power_of_two();

    let mut flattened_witness = vec![F::zero(); (num_registers_sup * num_rows_sup) as usize];

    worker.scope(num_rows, |scope, chunk| {
        // we take `chunks` registers and grab handle on the corresponding 
        for (i, f) in flattened_witness.chunks_mut(chunk * (num_registers_sup as usize)).enumerate() {
            let witness_iter = witness.iter();
            scope.spawn(move |_| {
                let start = chunk * i;
                let mut end = chunk * (i+1);
                if end > num_rows {
                    end = num_rows;
                }
                let range = end - start;
                for (reg_index, w) in witness_iter.clone().enumerate() {
                    for j in 0..range {
                        f[j*(num_registers_sup as usize) + reg_index] = w[j + start];
                    }
                }
            });
        }
    });
    
    let poly = Polynomial::<F, _>::from_values(flattened_witness)?;
    let witness_poly = poly.ifft(&worker);

    Ok(witness_poly)
}

impl<F: PrimeField> ARP<F> for SingleWitnessARPInstance<F> {
    fn from_instance(
        instance: IntoARPInstance<F>, 
        worker: &Worker
    ) -> Result<Self, SynthesisError> {
        let num_rows = instance.num_rows;
        let num_registers = instance.num_registers;

        let mut witness_poly = None;
        if let Some(w) = instance.witness {
            witness_poly = Some(make_single_witness_polymonial(w, &worker)?);
        }

        Ok(Self {
            witness_poly: witness_poly,
            num_rows: num_rows,
            num_registers: num_registers,
            constraints: instance.constraints,
            boundary_constraints: instance.boundary_constraints
        })
    }
}

impl<F: PrimeField> SingleWitnessARPInstance<F> {
    /// - make interpolating polynomial f
    /// - add masking coefficients for constraints
    /// - keep boundary constraints as it is
    pub fn route(&mut self) -> Result<(), SynthesisError> {
        let num_registers = self.num_registers as u64;
        let num_rows = self.num_rows as u64;

        let num_registers_sup = num_registers.next_power_of_two();
        let num_rows_sup = num_rows.next_power_of_two();

        let witness_domain = Domain::<F>::new_for_size(num_registers_sup * num_rows_sup)?;
        let rows_domain = Domain::<F>::new_for_size(num_rows_sup)?;

        fn remap_univariate_term<F: PrimeField>(
            term: &mut UnivariateTerm<F>,
            row_domain: &Domain<F>,
            column_domain: &Domain<F>,
        ) {
            let register_delta = match term.register {
                Register::Register(num) => {
                    num as u64
                },
                _ => {
                    unreachable!("Registers are now indistinguishable");
                }
            };

            let step_delta = match term.steps_difference {
                StepDifference::Steps(num) => {
                    num as u64
                },
                _ => {
                    unreachable!("Step differences are not masks yet");
                }
            };

            let mut mask = column_domain.generator.pow([step_delta]);
            let per_row_contrib = row_domain.generator.pow([register_delta]);

            mask.mul_assign(&per_row_contrib);

            term.steps_difference = StepDifference::Mask(mask);
        }

        fn remap_term<F: PrimeField>(
            term: &mut ConstraintTerm<F>,
            row_domain: &Domain<F>,
            column_domain: &Domain<F>,
        ) {
            match term {
                ConstraintTerm::Univariate(ref mut t) => {
                    remap_univariate_term(
                        t, 
                        &row_domain,
                        &column_domain
                    );
                },
                ConstraintTerm::Polyvariate(ref mut poly_term) => {
                    for mut t in poly_term.terms.iter_mut() {
                        remap_univariate_term(
                            &mut t, 
                            &row_domain,
                            &column_domain
                        );
                    }
                }
            }
        }

        fn remap_constraint<F: PrimeField>(
            constraint: &mut Constraint<F>,
            row_domain: &Domain<F>,
            column_domain: &Domain<F>,
        ) {
            for mut t in constraint.terms.iter_mut() {
                remap_term(
                    &mut t, 
                    &row_domain,
                    &column_domain
                );
            }
        }

        for mut c in self.constraints.iter_mut() {
            remap_constraint(
                &mut c, 
                &witness_domain,
                &rows_domain
            );
        }

        Ok(())
    }
}

#[test]
fn test_fib_conversion_into_arp() {
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
    let into_arp = test_tracer.into_arp();
    let worker = Worker::new();
    let mut arp = SingleWitnessARPInstance::<Fr>::from_instance(into_arp, &worker).expect("must work");
    arp.route().expect("must work");
    println!("Witness poly = {:?}", arp.witness_poly);

    for c in arp.constraints.iter() {
        println!("{:?}", c);
    }
}

