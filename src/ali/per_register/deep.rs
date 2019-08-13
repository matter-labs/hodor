use ff::PrimeField;

use crate::polynomials::*;
use crate::arp::*;
use crate::air::*;
use crate::fft::multicore::Worker;
use crate::SynthesisError;
use crate::domains::*;
use crate::precomputations::*;
use crate::transcript::Transcript;
use super::*;

impl<F: PrimeField> ALIInstance<F, PerRegisterARP> {
    pub fn calculate_deep<T: Transcript<F>>(
        &self,
        f_polys: &Vec<Polynomial<F, Coefficients>>,
        f_ldes: &Vec<Polynomial<F, Values>>,
        g_poly: &Polynomial<F, Coefficients>,
        g_lde: &Polynomial<F, Values>,
        transcript: &mut T,
        worker: &Worker
    ) -> Result<(Polynomial<F, Values>, Polynomial<F, Values>, Vec<F>, F), SynthesisError> {
        let z = transcript.get_challenge();

        let f_lde_size = (&f_ldes[0]).size();
        let g_lde_size = g_lde.size();
        let mut divisors_for_masks: IndexMap<StepDifference<F>, Polynomial<F, Values>> = IndexMap::new();

        let mut h1_lde = Polynomial::<F, Values>::new_for_size(f_lde_size)?;

        let mut f_at_z_m = vec![];

        println!("All masks length = {}", self.all_masks.len());

        for m in self.all_masks.iter() {
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

            let witness_value_at_z = f_polys[reg_num].evaluate_at(&worker, root);

            f_at_z_m.push(witness_value_at_z);
        
            let divisor = if let Some(d) = divisors_for_masks.get(&m.steps_difference) {
                d
            } else {
                let mut q_poly = Polynomial::<F, Coefficients>::new_for_size(2)?;
                q_poly.as_mut()[1] = F::one();
                q_poly.as_mut()[0].sub_assign(&root);
                let mut inverse_q_poly_coset_values = q_poly.evaluate_at_domain_for_degree_one(
                    &worker, 
                    f_lde_size as u64
                )?;

                inverse_q_poly_coset_values.batch_inversion(&worker)?;
                divisors_for_masks.insert(m.steps_difference, inverse_q_poly_coset_values);

                divisors_for_masks.get(&m.steps_difference).expect("is some now")
            };

            let mut f_minus_f_at_z = f_ldes[reg_num].clone();
            let mut minus_witness_value_at_z = witness_value_at_z;
            minus_witness_value_at_z.negate();
            f_minus_f_at_z.add_constant(&worker, &minus_witness_value_at_z);

            let alpha = transcript.get_challenge();

            f_minus_f_at_z.scale(&worker, alpha);
            f_minus_f_at_z.mul_assign(&worker, divisor);

            h1_lde.add_assign(&worker, &f_minus_f_at_z);

            // let alpha = transcript.get_challenge();
            // let beta = transcript.get_challenge();
            // let adjustment = self.max_constraint_power - 1;

            // let adj_poly = if adjustment == 0 {
            //     None
            // } else {
            //     let adj_poly = calculate_adjustment_polynomial_in_coset(
            //         adjustment,
            //         alpha,
            //         beta,
            //         &self.constraints_domain,
            //         &self.precomputations,
            //         &worker
            //     );

            //     Some(adj_poly)
            // };

            // let reg_num = match b_c.register {
            //     Register::Register(reg_number) => {
            //         reg_number
            //     },
            //     _ => {
            //         unreachable!();
            //     }
            // };
            // let mut witness_poly = (&witness[reg_num]).clone();
            // witness_poly.as_mut()[0].sub_assign(&b_c.value.expect("is some"));
            // let mut constraint_values = witness_poly.coset_lde(&worker, boundary_lde_factor as usize)?;
            // if let Some(adj) = adj_poly {
            //     constraint_values.mul_assign(&worker, &adj);
            // } else {
            //     // just apply alpha
            //     constraint_values.scale(&worker, alpha);
            // }

            // let divisors = self.boundary_constraint_divisors.get(&(b_c.at_row as u64)).expect("is some");
            // constraint_values.mul_assign(&worker, divisors);

            // g_values.add_assign(&worker, &constraint_values);
        }

        let mut q_poly = Polynomial::<F, Coefficients>::new_for_size(2)?;
        q_poly.as_mut()[1] = F::one();
        q_poly.as_mut()[0].sub_assign(&z);
        let mut inverse_q_poly_coset_values = q_poly.evaluate_at_domain_for_degree_one(
            &worker, 
            g_lde_size as u64
        )?;

        inverse_q_poly_coset_values.batch_inversion(&worker)?;

        let g_at_z = g_poly.evaluate_at(&worker, z);

        let mut h2_lde = g_lde.clone();
        let mut minus_g_at_z = g_at_z;
        minus_g_at_z.negate();

        h2_lde.add_constant(&worker, &minus_g_at_z);
        h2_lde.mul_assign(&worker, &inverse_q_poly_coset_values);

        Ok((h1_lde, h2_lde, f_at_z_m, g_at_z))
    }
}