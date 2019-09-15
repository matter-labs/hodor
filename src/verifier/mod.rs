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
use crate::fri::*;
use byteorder::{BigEndian, ByteOrder};

use indexmap::IndexSet as IndexSet;
use indexmap::IndexMap as IndexMap;
// use std::collections::{IndexSet, IndexMap};

/*

This module contains a stand-alone verifier. Verifier is composed from the instance definition and makes some precomputations
when initialized and all further verification calls are ammortized. Please mind that Stark verifier knows all the constraints,
so it requires linear memory to keep the constraints definition and linear time for evaluation of numerators for ALI step,
and at least logarithmic time to evaluation denominators in case of trivial "dense" constraints

Test in prover module containts a use example.

*/

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
        SynthesisError::Unsatisfied(format!("expecting value for term {:?} to be for register {} and step difference {:?}", univariate_term, reg_num, univariate_term.steps_difference))
    )?;

    let mut value = f_at_z_m.pow([univariate_term.power]);
    value.mul_assign(&univariate_term.coeff);

    Ok(value)
}
// ---------------------

pub struct InstanceProof<F: PrimeField, T: Transcript<F>, I: IOP<F>, P: FriProofPrototype<F, I>, PR: FriProof<F, I>, FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR>, A: ARPType> {
    pub f_at_z_m: Vec<F>, 
    pub f_iop_roots: Vec< < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput >,
    pub g_iop_root: < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput,

    pub f_queries: Vec<I::Query>,
    pub g_query: I::Query,

    pub h1_iop_roots: Vec< < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput >,
    pub h2_iop_roots: Vec< < <I::Tree as IopTree<F> >::Hasher as IopTreeHasher<F>>::HashOutput >,

    pub fri_proof_h1: PR,
    pub fri_proof_h2: PR,


    pub _marker_a: std::marker::PhantomData<A>,
    pub _marker_t: std::marker::PhantomData<T>,
    pub _marker_p: std::marker::PhantomData<P>,
    pub _marker_fri: std::marker::PhantomData<FRI>,
}

struct InstanceProofScratchSpace<F: PrimeField, T: Transcript<F>, I: IOP<F>, A: ARPType> {
    transcript: T,
    constraint_challenges: Vec<(F, F)>,
    boundary_constraint_challenges: Vec<(F, F)>,
    h1_poly_challenges: Vec<F>,

    _marker_a: std::marker::PhantomData<A>,
    _marker_i: std::marker::PhantomData<I>,
}

impl<F: PrimeField, T: Transcript<F>, I: IOP<F>, A: ARPType> InstanceProofScratchSpace<F, T, I, A> {
    fn new() -> Self {
        Self {
            transcript: T::new(),
            constraint_challenges: vec![],
            boundary_constraint_challenges: vec![],
            h1_poly_challenges: vec![],

            _marker_a: std::marker::PhantomData,
            _marker_i: std::marker::PhantomData,
        }
    }
}

pub struct Verifier<F: PrimeField, T: Transcript<F>, I: IOP<F>, P: FriProofPrototype<F, I>, PR: FriProof<F, I>, FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR>, A: ARPType> {
    instance: InstanceProperties<F>,
    max_constraint_power: u64,
    column_domain: Domain::<F>,
    constraints_domain: Domain::<F>,
    all_masks: IndexSet::<MaskProperties<F>>,
    all_boundary_constrained_registers: IndexSet::<Register>,
    constraints_batched_by_density: IndexMap::< ConstraintDensity, Vec<Constraint<F>> >,
    lde_factor: usize,
        
    _marker_a: std::marker::PhantomData<A>,
    _marker_i: std::marker::PhantomData<I>,
    _marker_t: std::marker::PhantomData<T>,
    _marker_p: std::marker::PhantomData<P>,
    _marker_fri: std::marker::PhantomData<FRI>,
}

impl<F: PrimeField, T: Transcript<F>, I: IOP<F>, P: FriProofPrototype<F, I>, PR: FriProof<F, I>, FRI: FriIop<F, IopType = I, ProofPrototype = P, Proof = PR>> Verifier<F, T, I, P, PR, FRI, PerRegisterARP> {
    pub fn new(
        instance: InstanceProperties<F>, 
        lde_factor: usize,
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

        for (_density, batch) in constraints_batched_by_density.iter() {
            for c in batch.iter()   {
                let constraint_power = c.degree;
                assert!(max_constraint_power >= constraint_power);
            }
        }

        Ok(Self {
            instance: instance,
            max_constraint_power: max_constraint_power,
            column_domain: column_domain,
            constraints_domain: constraints_domain,
            all_masks: all_masks,
            all_boundary_constrained_registers: all_boundary_constrained_registers,
            constraints_batched_by_density: constraints_batched_by_density,
            lde_factor,

            _marker_a: std::marker::PhantomData,
            _marker_i: std::marker::PhantomData,
            _marker_t: std::marker::PhantomData,
            _marker_p: std::marker::PhantomData,
            _marker_fri: std::marker::PhantomData,
        })
    }

    pub fn bytes_to_challenge_index<S: AsRef<[u8]>>(bytes: S, lde_size: usize, lde_factor: usize) -> usize {
        let as_ref = bytes.as_ref();
        let natural_x_index = BigEndian::read_u64(&as_ref[(as_ref.len() - 8)..]);

        let natural_x_index = natural_x_index as usize;
        let mut natural_x_index = natural_x_index % lde_size;
        if natural_x_index % lde_factor == 0 {
            natural_x_index += 1;
            natural_x_index = natural_x_index % lde_size;
        }

        if natural_x_index % 2 == 0 {
            natural_x_index += 1;
            natural_x_index = natural_x_index % lde_size;
        }

        natural_x_index
    }

    pub fn verify(
        &self,
        proof: &InstanceProof<F, T, I, P, PR, FRI, PerRegisterARP>
    ) -> Result<bool, SynthesisError> {
        let mut scratch_space: InstanceProofScratchSpace<F, T, I, PerRegisterARP> = InstanceProofScratchSpace::new();

        for r in proof.f_iop_roots.iter() {
            scratch_space.transcript.commit_bytes(r.as_ref());
        }

        for (_density, batch) in self.constraints_batched_by_density.iter() {
            for _c in batch.iter()   {
                let alpha = scratch_space.transcript.get_challenge();
                let beta = scratch_space.transcript.get_challenge();

                scratch_space.constraint_challenges.push((alpha, beta));
            }
        }

        for _ in self.instance.boundary_constraints.iter() {
            let alpha = scratch_space.transcript.get_challenge();
            let beta = scratch_space.transcript.get_challenge();

            scratch_space.boundary_constraint_challenges.push((alpha, beta));
        }

        // We've collected all challenges for G calculation and now can commit G
        scratch_space.transcript.commit_bytes(proof.g_iop_root.as_ref());

        let z = scratch_space.transcript.get_challenge();

        for _ in self.all_masks.iter() {
            let alpha = scratch_space.transcript.get_challenge();

            scratch_space.h1_poly_challenges.push(alpha);
        }

        // println!("Final root for h1 in verifier = {:?}", proof.h1_iop_roots.last().expect("there is one").as_ref());
        // println!("Final root for h1 in verifier = {:?}", proof.h2_iop_roots.last().expect("there is one").as_ref());

        scratch_space.transcript.commit_bytes(proof.h1_iop_roots.last().expect("there is one").as_ref());
        for el in proof.fri_proof_h1.get_final_coefficients().iter() {
            scratch_space.transcript.commit_field_element(&el);
        }

        scratch_space.transcript.commit_bytes(proof.h2_iop_roots.last().expect("there is one").as_ref());
        for el in proof.fri_proof_h2.get_final_coefficients().iter() {
            scratch_space.transcript.commit_field_element(&el);
        }

        let f_lde_size = self.column_domain.size * (self.lde_factor as u64);
        let g_lde_size = self.constraints_domain.size * (self.lde_factor as u64);

        let f_lde_domain = Domain::<F>::new_for_size(f_lde_size)?;
        let g_lde_domain = Domain::<F>::new_for_size(g_lde_size)?;

        println!("Getting challenge indexes");

        let x_challenge_index_h1 = Self::bytes_to_challenge_index(&scratch_space.transcript.get_challenge_bytes(), f_lde_size as usize, self.lde_factor);
        let x_challenge_index_h2 = Self::bytes_to_challenge_index(&scratch_space.transcript.get_challenge_bytes(), g_lde_size as usize, self.lde_factor);

        let mut f_ldes_at_x = vec![];

        if proof.f_queries.len() != self.instance.num_registers {
            return Err(SynthesisError::Unsatisfied(format!("expected number of registers = {}, proof contains queries only to {}", proof.f_queries.len(), self.instance.num_registers)));
        }

        if proof.f_queries.len() != proof.f_iop_roots.len() {
            return Err(SynthesisError::Unsatisfied(format!("received {} queries to witness oracles, but only {} roots", proof.f_queries.len(), proof.f_iop_roots.len())));
        }

        for (query, root) in proof.f_queries.iter().zip(proof.f_iop_roots.iter()) {
            if !I::verify_query(query, root) {
                return Ok(false);
            }
            if query.natural_index() != x_challenge_index_h1 {
                return Ok(false);
            }
            f_ldes_at_x.push(query.value());
        }

        println!("Simularing H1 part");

        let h_1_at_x = self.simulate_h1_from_f_at_z(
            &scratch_space.h1_poly_challenges, 
            x_challenge_index_h1, 
            &f_lde_domain, 
            &f_ldes_at_x, 
            &proof.f_at_z_m,
            z
        )?;

        println!("Calculating DEEP part");

        let g_at_z_from_verifier = self.calculate_g_at_z_from_f_at_z(
            &scratch_space,
            &proof,
            z
        )?;

        if !I::verify_query(&proof.g_query, &proof.g_iop_root) {
            return Ok(false);
        }
        if proof.g_query.natural_index() != x_challenge_index_h2 {
            return Ok(false);
        }

        let g_lde_at_x = proof.g_query.value();

        println!("Simularing H2 part");

        let h_2_at_x = self.simulate_h2_from_g_at_z(
            x_challenge_index_h2, 
            &g_lde_domain, 
            g_lde_at_x,
            z,
            g_at_z_from_verifier
        )?;

        // Now we need to check that H1 and H2 are indeed low degree polynomials
        let valid = FRI::verify_proof(
            &proof.fri_proof_h1,
            x_challenge_index_h1,
            h_1_at_x
        )?;

        if !valid {
            return Ok(false);
        }

        let valid = FRI::verify_proof(
            &proof.fri_proof_h2,
            x_challenge_index_h2,
            h_2_at_x
        )?;

        Ok(valid)
    }


    fn simulate_h1_from_f_at_z (
        &self,
        mask_challenges: &[F],
        natural_x_index: usize,
        f_lde_domain: &Domain<F>,
        f_ldes_at_x: &[F],
        f_at_z_m: &[F],
        z: F
    ) -> Result<F, SynthesisError> {
        let x = f_lde_domain.generator.pow(&[natural_x_index as u64]);

        let mut h_at_x = F::zero();

        // (f_i(x) - f(M*z)) / (x - M*z) = h_1(x)

        for ((m, f_at_z), alpha) in self.all_masks.iter()
                    .zip(f_at_z_m.iter())
                    .zip(mask_challenges.iter()) {
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
        scratch_space: &InstanceProofScratchSpace<F, T, I, PerRegisterARP>,
        proof: &InstanceProof<F, T, I, P, PR, FRI, PerRegisterARP>,
        z: F
    ) -> Result<F, SynthesisError> {
        let mut g_at_z = F::zero();
        let num_registers = self.instance.num_registers;

        let mut register_values_under_masks: Vec<IndexMap<StepDifference<F>, F>> = vec![IndexMap::new(); num_registers];

        for (m, f_at_z) in self.all_masks.iter()
                    .zip(proof.f_at_z_m.iter()) {
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

        let mut constraint_challenges_iter = scratch_space.constraint_challenges.iter();

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

        let mut boundary_constraint_challenges_iter = scratch_space.boundary_constraint_challenges.iter();

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