use ff::PrimeField;

use crate::arp::{ARPType, InstanceProperties, PerRegisterARP};
use crate::transcript::Transcript;
use crate::iop::*;
use crate::polynomials::*;
use crate::domains::*;
use crate::SynthesisError;
use crate::air::*;
use crate::ali::*;
use crate::arp::mappings::*;
use byteorder::{BigEndian, ByteOrder};

use indexmap::IndexSet as IndexSet;
use indexmap::IndexMap as IndexMap;
// use std::collections::{IndexSet, IndexMap};

pub struct Verifier<'a, F: PrimeField, T: Transcript<F>, I: IOP<'a, F>, A: ARPType> {
    instance: InstanceProperties<F>,
    max_constraint_power: u64,
    column_domain: Domain::<F>,
    constraints_domain: Domain::<F>,
    all_masks: IndexSet::<MaskProperties<F>>,
    f_at_z_m: Vec<F>, 
    all_boundary_constrained_registers: IndexSet::<Register>,
    constraints_batched_by_density: IndexMap::< ConstraintDensity, Vec<Constraint<F>> >,
    transcript: T,
    f_iop_roots: Vec< < <I::Tree as IopTree<'a, F> >::Hasher as IopTreeHasher<F>>::HashOutput >,
    g_iop_root: < <I::Tree as IopTree<'a, F> >::Hasher as IopTreeHasher<F>>::HashOutput,
    constraint_challenges: Vec<(F, F)>,
    boundary_constraint_challenges: Vec<(F, F)>,

    _marker: std::marker::PhantomData<A>
}

impl<'a, F: PrimeField, T: Transcript<F>, I: IOP<'a, F>> Verifier<'a, F, T, I, PerRegisterARP> {
    pub fn new(
        instance: InstanceProperties<F>, 
        f_at_z_m: Vec<F>,
        f_roots: Vec< < <I::Tree as IopTree<'a, F> >::Hasher as IopTreeHasher<F> >::HashOutput >,
        g_root: < <I::Tree as IopTree<'a, F> >::Hasher as IopTreeHasher<F> >::HashOutput
        ) -> Result<Self, SynthesisError> {

        let num_rows = instance.num_rows as u64;
        let num_rows_sup = num_rows.next_power_of_two();
        let column_domain = Domain::<F>::new_for_size(num_rows_sup)?;

        let mut instance = instance;

        for mut c in instance.constraints.iter_mut() {
            remap_constraint(
                &mut c, 
                &column_domain
            );
        }

        let mut transcript = T::new();
        for r in f_roots.iter() {
            transcript.commit_bytes(r.as_ref());
        }

        let mut all_masks: IndexSet::<MaskProperties<F>> = IndexSet::new();

        let mut max_constraint_power = 0u64;

        for c in instance.constraints.iter() {
            get_masks_from_constraint(&mut all_masks, c);
            if c.degree > max_constraint_power {
                max_constraint_power = c.degree;
            }
        }

        let constraint_power = max_constraint_power.next_power_of_two();

        let constraints_domain = Domain::<F>::new_for_size(constraint_power * num_rows_sup)?;

        let mut constraints_batched_by_density: IndexMap::< ConstraintDensity, Vec<Constraint<F>> > = IndexMap::new();

        for constraint in instance.constraints.iter() {
            if let Some(batch) = constraints_batched_by_density.get_mut(&constraint.density) {
                batch.push(constraint.clone());
            } else {
                constraints_batched_by_density.insert(constraint.density.clone(), vec![constraint.clone()]);
            }
        }

        let mut all_boundary_constrained_registers: IndexSet::<Register> = IndexSet::new();

        for b_c in instance.boundary_constraints.iter() {
            get_mask_from_boundary_constraint(&mut all_masks, b_c);

            if all_boundary_constrained_registers.get(&b_c.register).is_none() {
                all_boundary_constrained_registers.insert(b_c.register);
            }
        }

        // let mut all_boundary_constraint_steps: IndexSet<usize> = IndexSet::new();

        // for b_c in instance.boundary_constraints.iter() {
        //     if all_boundary_constraint_steps.get(&b_c.at_row).is_none() {
        //         all_boundary_constraint_steps.insert(b_c.at_row);
        //     }
        // }

        let mut constraint_challenges = vec![];

        for (_density, batch) in constraints_batched_by_density.iter() {
            for c in batch.iter()   {
                let constraint_power = c.degree;
                assert!(max_constraint_power >= constraint_power);
                let alpha = transcript.get_challenge();
                let beta = transcript.get_challenge();

                constraint_challenges.push((alpha, beta));
            }
        }

        let mut boundary_constraint_challenges = vec![];
        for _ in instance.boundary_constraints.iter() {
            let alpha = transcript.get_challenge();
            let beta = transcript.get_challenge();

            boundary_constraint_challenges.push((alpha, beta));
        }

        // We've collected all challenges for G calculation and now can commit G
        transcript.commit_bytes(g_root.as_ref());

        Ok(Self {
            instance: instance,
            max_constraint_power: max_constraint_power,
            column_domain: column_domain,
            constraints_domain: constraints_domain,
            all_masks: all_masks,
            f_at_z_m: f_at_z_m,
            all_boundary_constrained_registers: all_boundary_constrained_registers,
            constraints_batched_by_density: constraints_batched_by_density,
            transcript: transcript,
            f_iop_roots: f_roots,
            g_iop_root: g_root,
            constraint_challenges: constraint_challenges,
            boundary_constraint_challenges: boundary_constraint_challenges,

            _marker: std::marker::PhantomData
        })
    }

    fn simulate_h1_from_f_at_z (
        &self,
        transcript: T,
        natural_x_index: usize,
        f_lde_domain: &Domain<F>,
        f_ldes_at_x: &[F],
        z: F
    ) -> Result<F, SynthesisError> {
        let x = f_lde_domain.generator.pow(&[natural_x_index as u64]);

        let mut h_at_x = F::zero();

        let mut transcript = transcript;

        // (f_i(x) - f(M*z)) / (x - M*z) = h_1(x)

        for (m, f_at_z) in self.all_masks.iter()
                    .zip(self.f_at_z_m.iter()) {
            let mut root = match m.steps_difference {
                StepDifference::Mask(mask) => {
                    mask
                },
                _ => {
                    unreachable!();
                }
            };
            root.mul_assign(&z);

            let reg_num = match m.register {
                Register::Register(reg_number) => {
                    reg_number
                },
                _ => {
                    unreachable!();
                }
            };

            let f_at_x = f_ldes_at_x[reg_num];

            let mut num = f_at_x;
            num.sub_assign(&f_at_z);

            let mut den = x;
            den.sub_assign(&root);

            let den_inv = den.inverse().ok_or(
                SynthesisError::DivisionByZero(format!("no inverse for mask {:?} at z = {}", m, z))
            )?;

            num.mul_assign(&den_inv);

            let alpha = transcript.get_challenge();
            num.mul_assign(&alpha);

            h_at_x.add_assign(&num);
        }

        Ok(h_at_x)
    }

    fn simulate_h2_from_g_at_z (
        &self,
        natural_x_index: usize,
        g_lde_domain: &Domain<F>,
        g_lde_at_x: F,
        z: F,
        g_at_z: F
    ) -> Result<F, SynthesisError> {
        let x = g_lde_domain.generator.pow(&[natural_x_index as u64]);

        // (g(x) - g(z)) / (x - z) = h_2(x)

        let mut h_at_x = g_lde_at_x;
        h_at_x.sub_assign(&g_at_z);

        let mut den = x;
        den.sub_assign(&z);

        let den_inv = den.inverse().ok_or(
            SynthesisError::DivisionByZero(format!("no inverse at z = {}", z))
        )?;

        h_at_x.mul_assign(&den_inv);

        Ok(h_at_x)
    }

    fn calculate_g_at_z_from_f_at_z (
        &self,
        z: F
    ) -> Result<F, SynthesisError> {
        let mut g_at_z = F::zero();
        let num_registers = self.instance.num_registers;

        // ---------------------
        fn evaluate_constraint_on_f_at_z_m<F: PrimeField>(
            constraint: &Constraint<F>,
            witness: &Vec<IndexMap<StepDifference<F>, F>>,
        ) -> Result<F, SynthesisError> {
            let mut value = constraint.constant_term;
            for term in constraint.terms.iter() {
                let v = evaluate_term_on_f_at_z_m(&term, witness)?;
                value.add_assign(&v);
            }
            
            Ok(value)
        }

        fn evaluate_term_on_f_at_z_m<F: PrimeField>(
            term: &ConstraintTerm<F>,
            witness: &Vec<IndexMap<StepDifference<F>, F>>,
        ) -> Result<F, SynthesisError> {
            match term {
                ConstraintTerm::Univariate(ref t) => {
                    evaluate_univariate_term_on_f_at_z_m(
                        t, 
                        &witness,
                    )
                },
                ConstraintTerm::Polyvariate(ref poly_term) => {
                    let mut result = F::one();
                    for t in poly_term.terms.iter() {
                        let v = evaluate_univariate_term_on_f_at_z_m(
                            &t, 
                            &witness,
                        )?;
                        result.mul_assign(&v);
                    }

                    result.mul_assign(&poly_term.coeff);

                    Ok(result)
                }
            }
        }

        fn evaluate_univariate_term_on_f_at_z_m<F: PrimeField>(
            univariate_term: &UnivariateTerm<F>,
            witness: &Vec<IndexMap<StepDifference<F>, F>>,
        ) -> Result<F, SynthesisError> {
            // let mask = match univariate_term.steps_difference {
            //     StepDifference::Mask(mask) => mask,
            //     _ => unreachable!()
            // };

            let reg_num = match univariate_term.register {
                Register::Register(reg_num) => reg_num,
                _ => unreachable!()
            };

            let f_at_z_m = witness[reg_num].get(&univariate_term.steps_difference).ok_or(
                SynthesisError::Unsatisfied(format!("Expecting value for term {:?} to be for register {} and step difference {:?}", univariate_term, reg_num, univariate_term.steps_difference))
            )?;

            let mut value = f_at_z_m.pow([univariate_term.power]);
            value.mul_assign(&univariate_term.coeff);

            Ok(value)
        }
        // ---------------------

        let mut register_values_under_masks: Vec<IndexMap<StepDifference<F>, F>> = vec![IndexMap::new(); num_registers];

        for (m, f_at_z) in self.all_masks.iter()
                    .zip(self.f_at_z_m.iter()) {
            let mut root = match m.steps_difference {
                StepDifference::Mask(mask) => {
                    mask
                },
                _ => {
                    unreachable!();
                }
            };
            root.mul_assign(&z);

            let reg_num = match m.register {
                Register::Register(reg_number) => {
                    reg_number
                },
                _ => {
                    unreachable!();
                }
            };

            register_values_under_masks[reg_num].insert(m.steps_difference, *f_at_z);
        }

        let mut constraint_challenges_iter = self.constraint_challenges.iter();

        for (density, batch) in self.constraints_batched_by_density.iter() {
            let (inverse_divisor, _) = match density {
                ConstraintDensity::Dense(dense) => {
                    inverse_divisor_for_dense_constraint(
                        z, 
                        &self.column_domain, 
                        &self.constraints_domain, 
                        dense.clone(),
                        self.instance.num_rows as u64
                )? 
                },
                _ => {
                    unimplemented!();
                }
            };
            
            for c in batch.iter() {
                let constraint_power = c.degree;
                assert!(self.max_constraint_power >= constraint_power);

                let (alpha, beta) = constraint_challenges_iter.next().expect("have enough challenges");

                let mut value_at_z = evaluate_constraint_on_f_at_z_m(
                    c, 
                    &register_values_under_masks
                )?;

                let adjustment = self.max_constraint_power - constraint_power;

                if adjustment == 0 {
                    value_at_z.mul_assign(&alpha);
                } else {
                    // adjustment = alpha * z^(adj) + beta

                    let mut adj = z.pow([adjustment]);
                    adj.mul_assign(&alpha);
                    adj.add_assign(&beta);
                    
                    value_at_z.mul_assign(&adj);
                };

                value_at_z.mul_assign(&inverse_divisor);

                g_at_z.add_assign(&value_at_z);
            }
        }

        let mut boundary_constraint_challenges_iter = self.boundary_constraint_challenges.iter();

        let current_step = StepDifference::Mask(F::one());

        for b_c in self.instance.boundary_constraints.iter() {
            let row = b_c.at_row;

            let adjustment = self.max_constraint_power - 1;

            let (alpha, beta) = boundary_constraint_challenges_iter.next().expect("have enough challenges");

            let reg_num = match b_c.register {
                Register::Register(reg_number) => {
                    reg_number
                },
                _ => {
                    unreachable!();
                }
            };

            let mut value_at_z = *register_values_under_masks[reg_num].get(&current_step).ok_or(
                SynthesisError::Unsatisfied(format!("Expecting value for boundary constraint {:?} to be for register {} and no step difference", b_c, reg_num))
            )?;

            let expected_value = b_c.value.ok_or(
                SynthesisError::Unsatisfied(format!("Expecting value for boundary constraint {:?} for register {} and at row {}", b_c, reg_num, row))
            )?;

            value_at_z.sub_assign(&expected_value);

            let root = self.column_domain.generator.pow([row as u64]);

            let mut divisor = z;
            divisor.sub_assign(&root);

            let inverse_divisor = divisor.inverse().ok_or(
                SynthesisError::DivisionByZero(format!("no inverse for boundary constraint {:?} at z = {}", b_c, z))
            )?;

            if adjustment == 0 {
                value_at_z.mul_assign(&alpha);
            } else {
                // adjustment = alpha * z^(adj) + beta

                let mut adj = z.pow([adjustment]);
                adj.mul_assign(&alpha);
                adj.add_assign(&beta);
                
                value_at_z.mul_assign(&adj);
            };

            value_at_z.mul_assign(&inverse_divisor);

            g_at_z.add_assign(&value_at_z);
        }

        Ok(g_at_z)
    }
}

// such calls most likely will have start at 0 and num_steps = domain_size - 1
fn inverse_divisor_for_dense_constraint<F: PrimeField> (
    x: F,
    column_domain: &Domain<F>,
    _evaluation_domain: &Domain<F>,
    dense_constraint: DenseConstraint,
    num_rows: u64,
) -> Result<(F, usize), SynthesisError> {
    let start_at = dense_constraint.start_at;
    let span = dense_constraint.span as u64;
    let mut divisor_degree = column_domain.size as usize;
    let divisor_domain_size = column_domain.size;
    divisor_degree -= start_at as usize;
    divisor_degree -= (divisor_domain_size - num_rows) as usize;
    divisor_degree -= span as usize;

    let mut q_at_x = x.pow(&[divisor_domain_size]);
    q_at_x.sub_assign(&F::one());

    let mut inverse_divisor = q_at_x.inverse().ok_or(
        SynthesisError::DivisionByZero(format!("no inverse for dense constraint {:?} at x = {}", dense_constraint, x))
    )?;

    let roots_generator = column_domain.generator;

    let mut root = F::one();
    for _ in 0..start_at {
        let mut tmp = x;
        tmp.sub_assign(&root);
        inverse_divisor.mul_assign(&tmp);
        root.mul_assign(&roots_generator);                
    }
    
    let last_step = num_rows - span;
    let mut root = roots_generator.pow([last_step]);
    for _ in last_step..divisor_domain_size {
        let mut tmp = x;
        tmp.sub_assign(&root);
        inverse_divisor.mul_assign(&tmp);
        root.mul_assign(&roots_generator);
    }

    Ok((inverse_divisor, divisor_degree))
}

#[test]
fn test_fib_full_verifier() {
    use ff::Field;
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::fft::multicore::Worker;
    use crate::transcript::Transcript;
    use crate::arp::*;
    use crate::ali::*;
    use crate::iop::blake2s_trivial_iop::TrivialBlake2sIOP;
    use crate::iop::blake2s_trivial_iop::Blake2sIopTree;
    use crate::fri::*;
    use crate::transcript::*;
    use crate::ali::per_register::*;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let lde_factor = 16;
    let mut transcript = Blake2sTranscript::new();

    let worker = Worker::new();

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let (witness, props) = test_tracer.into_arp();
    let witness = witness.expect("some witness");
    // println!("Witness = {:?}", witness);

    let is_satisfied = ARPInstance::<Fr, PerRegisterARP>::is_satisfied(&props, &witness, &worker);
    assert!(is_satisfied.is_ok());

    let arp = ARPInstance::<Fr, PerRegisterARP>::from_instance(props.clone(), &worker).expect("must work");

    let witness_polys = arp.calculate_witness_polys(witness, &worker).expect("must work");

    let f_ldes: Vec<_> = witness_polys.iter().map(|w| {
        w.clone().lde(&worker, lde_factor).expect("must work")
    }).collect();

    let f_oracles: Vec<_> = f_ldes.iter().map(|l|
        Blake2sIopTree::create(l.as_ref())
    ).collect(); 

    for o in f_oracles.iter() {
        transcript.commit_bytes(&o.get_root()[..]);
    }

    let ali = ALIInstance::from_arp(arp, &worker).expect("is some");

    let g_poly_interpolant = ali.calculate_g(&mut transcript, witness_polys.clone(), &worker).expect("is some");

    let g_lde = g_poly_interpolant.clone().lde(&worker, lde_factor).expect("is something");

    let g_oracle = Blake2sIopTree::create(g_lde.as_ref());
    transcript.commit_bytes(&g_oracle.get_root());

    let (h1_lde, h2_lde, f_at_z_m, g_at_z) = ali.calculate_deep(
        &witness_polys,
        &f_ldes,
        &g_poly_interpolant,
        &g_lde,
        &mut transcript,
        &worker
    ).expect("must work");

    // type FriProver = NaiveFriIop::<'_, Fr, TrivialBlake2sIOP<Fr>>;
    //  as FRIIOP<'_, Fr>;

    let fri_final_poly_degree = 1;

    let h1_fri_proof_proto = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde(&h1_lde, lde_factor, fri_final_poly_degree, &worker).expect("must work");
    let h2_fri_proof_proto = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::proof_from_lde(&h2_lde, lde_factor, fri_final_poly_degree, &worker).expect("must work");
    
    let challenge_root_h1 = h1_fri_proof_proto.final_root;

    let h1_lde_size = h1_lde.size();

    let natural_x_index = BigEndian::read_u64(&challenge_root_h1[(challenge_root_h1.len() - 8)..]);

    let natural_x_index = natural_x_index as usize;
    let mut natural_x_index_h1 = natural_x_index % h1_lde_size;
    if natural_x_index_h1 % lde_factor == 0 {
        natural_x_index_h1 += 1;
        natural_x_index_h1 = natural_x_index_h1 % h1_lde_size;
    }

    if natural_x_index_h1 % 2 == 0 {
        natural_x_index_h1 += 1;
        natural_x_index_h1 = natural_x_index_h1 % h1_lde_size;
    }

    let h1_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::prototype_into_proof(h1_fri_proof_proto, &h1_lde, natural_x_index_h1).expect("must work");

    let challenge_root_h2 = h2_fri_proof_proto.final_root;

    let h2_lde_size = h2_lde.size();

    let natural_x_index = BigEndian::read_u64(&challenge_root_h2[(challenge_root_h2.len() - 8)..]);

    let natural_x_index = natural_x_index as usize;
    let mut natural_x_index_h2 = natural_x_index % h2_lde_size;
    if natural_x_index_h2 % lde_factor == 0 {
        natural_x_index_h2 += 1;
        natural_x_index_h2 = natural_x_index_h2 % h2_lde_size;
    }

    if natural_x_index_h2 % 2 == 0 {
        natural_x_index_h2 += 1;
        natural_x_index_h2 = natural_x_index_h2 % h2_lde_size;
    }

    let h2_fri_proof = NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>>::prototype_into_proof(h2_fri_proof_proto, &h2_lde, natural_x_index_h2).expect("must work");

    // All prover work is complete here 

    let mut f_roots = vec![];

    for o in f_oracles.iter() {
        f_roots.push(o.get_root());
    }

    let g_root = g_oracle.get_root();

    let mut verifier = Verifier::<Fr, Blake2sTranscript<Fr>, TrivialBlake2sIOP<Fr>, PerRegisterARP>::new(
        props, 
        f_at_z_m,
        f_roots,
        g_root,
    ).expect("some verifier");

    let mut f_ldes_at_x = vec![];
    for f in f_ldes.iter() {
        f_ldes_at_x.push(f.as_ref()[natural_x_index_h1]);
    }

    let z = verifier.transcript.get_challenge();

    let f_lde_domain = Domain::<Fr>::new_for_size(f_ldes[0].size() as u64).expect("some domain");

    let h_1_at_x = verifier.simulate_h1_from_f_at_z(
        verifier.transcript.clone(), 
        natural_x_index_h1, 
        &f_lde_domain, 
        &f_ldes_at_x, 
        z
    ).expect("some h_1 value");

    assert_eq!(h_1_at_x, h1_lde.as_ref()[natural_x_index_h1], "h_1 simulation failed");

    let g_at_z_from_verifier = verifier.calculate_g_at_z_from_f_at_z(z).expect("some g at z");

    assert_eq!(g_at_z, g_at_z_from_verifier, "g at z is not the same in prover and verifier");

    let g_lde_domain = Domain::<Fr>::new_for_size(g_lde.size() as u64).expect("some domain");

    let g_lde_at_x = g_lde.as_ref()[natural_x_index_h2];

    let h_2_at_x = verifier.simulate_h2_from_g_at_z(
        natural_x_index_h2, 
        &g_lde_domain, 
        g_lde_at_x,
        z,
        g_at_z_from_verifier
    ).expect("some h_2 value");

    assert_eq!(h_2_at_x, h2_lde.as_ref()[natural_x_index_h2], "h_2 simulation failed");

    // Now we need to check that H1 and H2 are indeed low degree polynomials

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::verify_proof(
        &h1_fri_proof,
        natural_x_index_h1,
        h_1_at_x
    ).expect("must work");

    assert!(valid);

    let valid = <NaiveFriIop::<Fr, TrivialBlake2sIOP<Fr>> as FRIIOP<'_, Fr>>::verify_proof(
        &h2_fri_proof,
        natural_x_index_h2,
        h_2_at_x
    ).expect("must work");

    assert!(valid);
}