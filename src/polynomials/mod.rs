use ff::PrimeField;
// use crate::plonk::domains::*;
use crate::domains::*;
use crate::fft::multicore::*;
use crate::fft::*;
use crate::{SynthesisError};
// use crate::plonk::fft::with_precomputation;
// use crate::plonk::fft::with_precomputation::FftPrecomputations;
// use crate::plonk::fft::*;
// use crate::worker::*;

use crate::polynomials::cooley_tukey_ntt::*;
// use crate::plonk::fft::cooley_tukey_ntt;
// use crate::plonk::fft::cooley_tukey_ntt::partial_reduction;
// use crate::plonk::fft::cooley_tukey_ntt::CTPrecomputations;

// use crate::plonk::transparent_engine::PartialTwoBitReductionField;

pub trait PolynomialForm: Sized + Copy + Clone {}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Coefficients {}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Values {}

impl PolynomialForm for Coefficients {}
impl PolynomialForm for Values {}

// TODO: Enforce bitreversed values as a separate form

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Polynomial<F: PrimeField, P: PolynomialForm> {
    coeffs: Vec<F>,
    pub exp: u32,
    pub omega: F,
    pub omegainv: F,
    pub geninv: F,
    pub minv: F,
    _marker: std::marker::PhantomData<P>,
}

impl<F: PrimeField, P: PolynomialForm> Polynomial<F, P> {
    pub fn size(&self) -> usize {
        self.coeffs.len()
    }

    pub fn as_ref(&self) -> &[F] {
        &self.coeffs
    }

    pub fn as_mut(&mut self) -> &mut [F] {
        &mut self.coeffs
    }

    pub fn into_coeffs(self) -> Vec<F> {
        self.coeffs
    }

    pub fn distribute_powers(&mut self, worker: &Worker, g: F) {
        distribute_powers(&mut self.coeffs, &worker, g);
    }

    pub fn reuse_allocation<PP: PolynomialForm>(&mut self, other: &Polynomial<F, PP>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());
        self.coeffs.copy_from_slice(&other.coeffs);
    }

    pub fn bitreverse_enumeration(&mut self, worker: &Worker) {
        let total_len = self.coeffs.len();
        let log_n = self.exp as usize;
        if total_len <= worker.cpus {
            for j in 0..total_len {
                let rj = cooley_tukey_ntt::bitreverse(j, log_n);
                if j < rj {
                    self.coeffs.swap(j, rj);
                }
            }

            return;
        }

        let r = &mut self.coeffs[..] as *mut [F];

        let to_spawn = worker.cpus;

        let chunk = Worker::chunk_size_for_num_spawned_threads(total_len, to_spawn);

        // while it's unsafe we don't care cause swapping is always one to one

        worker.scope(0, |scope, _| {
            for thread_id in 0..to_spawn {
                let r = unsafe { &mut *r };
                scope.spawn(move |_| {
                    let start = thread_id * chunk;
                    let end = if start + chunk <= total_len {
                        start + chunk
                    } else {
                        total_len
                    };
                    for j in start..end {
                        let rj = cooley_tukey_ntt::bitreverse(j, log_n);
                        if j < rj {
                            r.swap(j, rj);
                        }
                    }
                });
            }
        });
    }

    pub fn scale(&mut self, worker: &Worker, g: F) {
        if g == F::one() {
            return;
        }

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        v.mul_assign(&g);
                    }
                });
            }
        });
    }

    pub fn negate(&mut self, worker: &Worker) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        v.negate();
                    }
                });
            }
        });
    }

    pub fn map<M: Fn(&mut F) -> () + Send + Copy>(&mut self, worker: &Worker, func: M) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        func(v);
                    }
                });
            }
        });
    }

    pub fn map_indexed<M: Fn(usize, &mut F) -> () + Send + Copy>(
        &mut self,
        worker: &Worker,
        func: M,
    ) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (chunk_idx, v) in self.coeffs.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let base = chunk_idx * chunk;
                    for (in_chunk_idx, v) in v.iter_mut().enumerate() {
                        func(base + in_chunk_idx, v);
                    }
                });
            }
        });
    }

    pub fn pad_by_factor(&mut self, factor: usize) -> Result<(), SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(());
        }
        let next_power_of_two = factor.next_power_of_two();
        if factor != next_power_of_two {
            return Err(SynthesisError::Unsatisfied("unsatisfiable constraint system".to_string()));
        }

        let new_size = self.coeffs.len() * factor;
        self.coeffs.resize(new_size, F::zero());

        let domain = Domain::new_for_size(new_size as u64)?;
        self.exp = domain.power_of_two as u32;
        let m = domain.size as usize;
        self.omega = domain.generator;
        self.omegainv = self.omega.inverse().unwrap();
        self.minv = F::from_str(&format!("{}", m)).unwrap().inverse().unwrap();

        Ok(())
    }

    pub fn pad_to_size(&mut self, new_size: usize) -> Result<(), SynthesisError> {
        if new_size < self.coeffs.len() {
            return Err(SynthesisError::PolynomialDegreeTooLarge);
        }
        let next_power_of_two = new_size.next_power_of_two();
        if new_size != next_power_of_two {
            return Err(SynthesisError::Unsatisfied("unsatisfiable constraint system".to_string()));
        }
        self.coeffs.resize(new_size, F::zero());

        let domain = Domain::new_for_size(new_size as u64)?;
        self.exp = domain.power_of_two as u32;
        let m = domain.size as usize;
        self.omega = domain.generator;
        self.omegainv = self.omega.inverse().unwrap();
        self.minv = F::from_str(&format!("{}", m)).unwrap().inverse().unwrap();

        Ok(())
    }

    pub fn pad_to_domain(&mut self) -> Result<(), SynthesisError> {
        let domain = Domain::<F>::new_for_size(self.coeffs.len() as u64)?;
        self.coeffs.resize(domain.size as usize, F::zero());

        Ok(())
    }

    pub fn clone_padded_to_domain(&self) -> Result<Self, SynthesisError> {
        let mut padded = self.clone();
        let domain = Domain::<F>::new_for_size(self.coeffs.len() as u64)?;
        padded.coeffs.resize(domain.size as usize, F::zero());

        Ok(padded)
    }

    pub fn trim_to_degree(&mut self, degree: usize) -> Result<(), SynthesisError> {
        let size = self.coeffs.len();
        if size <= degree + 1 {
            return Ok(());
        }
        self.coeffs.truncate(degree + 1);
        self.coeffs.resize(size, F::zero());

        Ok(())
    }
}

impl<F: PrimeField> Polynomial<F, Coefficients> {
    pub fn new_for_size(size: usize) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
        let coeffs = vec![F::zero(); size];

        Self::from_coeffs(coeffs)
    }

    pub fn from_coeffs(mut coeffs: Vec<F>) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
        let coeffs_len = coeffs.len();

        let domain = Domain::new_for_size(coeffs_len as u64)?;
        let exp = domain.power_of_two as u32;
        let m = domain.size as usize;
        let omega = domain.generator;

        coeffs.resize(m, F::zero());

        Ok(Polynomial::<F, Coefficients> {
            coeffs: coeffs,
            exp: exp,
            omega: omega,
            omegainv: omega.inverse().unwrap(),
            geninv: F::multiplicative_generator().inverse().unwrap(),
            minv: F::from_str(&format!("{}", m)).unwrap().inverse().unwrap(),
            _marker: std::marker::PhantomData,
        })
    }

    pub fn from_roots(
        roots: Vec<F>,
        worker: &Worker,
    ) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
        let coeffs_len = roots.len() + 1;

        let domain = Domain::<F>::new_for_size(coeffs_len as u64)?;
        let num_threads = worker.cpus;

        // vector of vectors of polynomial coefficients for subproducts
        let mut subterms = vec![vec![]; num_threads];

        worker.scope(roots.len(), |scope, chunk| {
            for (r, s) in roots.chunks(chunk).zip(subterms.chunks_mut(1)) {
                scope.spawn(move |_| {
                    for r in r.iter() {
                        if s[0].len() == 0 {
                            let mut tmp = *r;
                            tmp.negate();
                            s[0] = vec![tmp, F::one()];
                        } else {
                            let mut tmp = Vec::with_capacity(s[0].len() + 1);
                            tmp.push(F::zero());
                            tmp.extend(s[0].clone());
                            for (c, n) in s[0].iter().zip(tmp.iter_mut()) {
                                let mut t = *c;
                                t.mul_assign(&r);
                                n.sub_assign(&t);
                            }
                            s[0] = tmp;
                        }
                    }
                });
            }
        });

        // now we have subproducts in a coefficient form

        let mut result: Option<Polynomial<F, Values>> = None;
        let result_len = domain.size as usize;

        for s in subterms.into_iter() {
            if s.len() == 0 {
                continue;
            }
            let t = Polynomial::<F, Coefficients>::from_coeffs(s)?;
            let factor = result_len / t.size();
            let t = t.lde(&worker, factor)?;
            if let Some(res) = result.as_mut() {
                res.mul_assign(&worker, &t);
            } else {
                result = Some(t);
            }
        }

        let result = result.expect("is some");
        let result = result.ifft(&worker);

        Ok(result)
    }

    pub fn evaluate_at_domain_for_degree_one(
        &self,
        worker: &Worker,
        domain_size: u64,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        assert_eq!(self.coeffs.len(), 2);
        let alpha = self.coeffs[1];
        let c = self.coeffs[0];

        let domain = Domain::<F>::new_for_size(domain_size)?;

        let mut result = vec![F::zero(); domain.size as usize];
        let g = domain.generator;
        worker.scope(result.len(), |scope, chunk| {
            for (i, v) in result.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = g.pow(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        let mut tmp = alpha;
                        tmp.mul_assign(&u);
                        tmp.add_assign(&c);
                        *v = tmp;
                        u.mul_assign(&g);
                    }
                });
            }
        });

        Polynomial::from_values(result)
    }

    pub fn coset_evaluate_at_domain_for_degree_one(
        &self,
        worker: &Worker,
        domain_size: u64,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        assert_eq!(self.coeffs.len(), 2);
        let alpha = self.coeffs[1];
        let c = self.coeffs[0];

        let domain = Domain::<F>::new_for_size(domain_size)?;

        let mut result = vec![F::zero(); domain.size as usize];
        let g = domain.generator;
        worker.scope(result.len(), |scope, chunk| {
            for (i, v) in result.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = g.pow(&[(i * chunk) as u64]);
                    u.mul_assign(&F::multiplicative_generator());
                    for v in v.iter_mut() {
                        let mut tmp = alpha;
                        tmp.mul_assign(&u);
                        tmp.add_assign(&c);
                        *v = tmp;
                        u.mul_assign(&g);
                    }
                });
            }
        });

        Polynomial::from_values(result)
    }

    // pub fn sparse_distribute(&mut self, factor: usize, worker: &Worker) -> Result<(), SynthesisError> {
    //     if factor == 1 {
    //         return Ok(());
    //     }
    //     let next_power_of_two = factor.next_power_of_two();
    //     if factor != next_power_of_two {
    //         return Err(SynthesisError::Error);
    //     }

    //     let new_size = self.coeffs.len() * factor;
    //     let new_coeffs = vec![F::zero(); new_size];
    //     let old_coeffs = std::mem::replace(&mut self.coeffs, new_coeffs);

    //     // now we need to interleave the coefficients
    //     worker.scope(old_coeffs.len(), |scope, chunk| {
    //         for (old, new) in old_coeffs.chunks(chunk)
    //                         .zip(self.coeffs.chunks_mut(chunk*factor)) {
    //             scope.spawn(move |_| {
    //                 for (j, old_coeff) in old.iter().enumerate() {
    //                     new[j*factor] = *old_coeff;
    //                 }
    //             });
    //         }
    //     });

    //     let domain = Domain::new_for_size(new_size as u64)?;
    //     self.exp = domain.power_of_two as u32;
    //     let m = domain.size as usize;
    //     self.omega = domain.generator;
    //     self.omegainv = self.omega.inverse().unwrap();
    //     self.minv = F::from_str(&format!("{}", m)).unwrap().inverse().unwrap();

    //     Ok(())
    // }

    // pub fn extend(&mut self, factor: usize, _worker: &Worker) -> Result<(), SynthesisError> {
    //     if factor == 1 {
    //         return Ok(());
    //     }
    //     let next_power_of_two = factor.next_power_of_two();
    //     if factor != next_power_of_two {
    //         return Err(SynthesisError::Error);
    //     }

    //     let new_size = self.coeffs.len() * factor;
    //     self.coeffs.resize(new_size, F::zero());

    //     Ok(())
    // }

    #[inline(always)]
    pub fn break_into_multiples(
        self,
        size: usize,
    ) -> Result<Vec<Polynomial<F, Coefficients>>, SynthesisError> {
        let mut coeffs = self.coeffs;

        let (mut num_parts, last_size) = if coeffs.len() % size == 0 {
            let num_parts = coeffs.len() / size;

            (num_parts, 0)
        } else {
            let num_parts = coeffs.len() / size;
            let last_size = coeffs.len() - num_parts * size;

            (num_parts, last_size)
        };

        let mut rev_results = Vec::with_capacity(num_parts);

        if last_size != 0 {
            let drain = coeffs.split_off(coeffs.len() - last_size);
            rev_results.push(drain);
            num_parts -= 1;
        }

        for _ in 0..num_parts {
            let drain = coeffs.split_off(coeffs.len() - size);
            rev_results.push(drain);
        }

        let mut results = Vec::with_capacity(num_parts);

        for c in rev_results.into_iter().rev() {
            let poly = Polynomial::<F, Coefficients>::from_coeffs(c)?;
            results.push(poly);
        }

        Ok(results)
    }

    #[inline(always)]
    pub fn lde(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        self.lde_using_multiple_cosets(worker, factor)
        // self.filtering_lde(worker, factor)
    }

    #[inline(always)]
    pub fn coset_lde(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        self.coset_lde_using_multiple_cosets(worker, factor)
        // self.filtering_coset_lde(worker, factor)
    }

    pub fn filtering_lde(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(self.fft(&worker));
        }
        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::new_for_size(new_size as u64)?;

        let mut lde = self.coeffs;
        lde.resize(new_size as usize, F::zero());
        best_lde(
            &mut lde,
            worker,
            &domain.generator,
            domain.power_of_two as u32,
            factor,
        );

        Polynomial::from_values(lde)
    }

    pub fn lde_using_multiple_cosets_naive(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(self.fft(&worker));
        }
        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::new_for_size(new_size as u64)?;

        let mut results = vec![];

        let mut coset_generator = F::one();

        let one = F::one();

        for _ in 0..factor {
            let coeffs = self.clone();
            let lde = if coset_generator == one {
                coeffs.fft(&worker)
            } else {
                coeffs.coset_fft_for_generator(&worker, coset_generator)
            };

            results.push(lde.into_coeffs());
            coset_generator.mul_assign(&domain.generator);
        }

        let mut final_values = vec![F::zero(); new_size];

        let results_ref = &results;

        worker.scope(final_values.len(), |scope, chunk| {
            for (i, v) in final_values.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = i * chunk;
                    for v in v.iter_mut() {
                        let coset_idx = idx % factor;
                        let element_idx = idx / factor;
                        *v = results_ref[coset_idx][element_idx];

                        idx += 1;
                    }
                });
            }
        });

        Polynomial::from_values(final_values)
    }

    pub fn lde_using_multiple_cosets(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(self.fft(&worker));
        }

        let num_cpus = worker.cpus;
        let num_cpus_hint = if num_cpus <= factor {
            Some(1)
        } else {
            let threads_per_coset = factor / num_cpus;
            // TODO: it's not immediately clean if having more threads than (virtual) cores benefits
            // over underutilization of some (virtual) cores
            // let mut threads_per_coset = factor / num_cpus;
            // if factor % num_cpus != 0 {
            //     if (threads_per_coset + 1).is_power_of_two() {
            //         threads_per_coset += 1;
            //     }
            // }
            Some(threads_per_coset)
        };

        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::<F>::new_for_size(new_size as u64)?;

        let mut results = vec![self.coeffs; factor];

        let coset_omega = domain.generator;
        let this_domain_omega = self.omega;

        let log_n = self.exp;

        worker.scope(results.len(), |scope, chunk| {
            for (i, r) in results.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut coset_generator = coset_omega.pow(&[i as u64]);
                    for r in r.iter_mut() {
                        if coset_generator != F::one() {
                            distribute_powers_serial(&mut r[..], coset_generator);
                            // distribute_powers(&mut r[..], &worker, coset_generator);
                        }
                        best_fft(
                            &mut r[..],
                            &worker,
                            &this_domain_omega,
                            log_n,
                            // num_cpus_hint,
                            None,
                        );
                        coset_generator.mul_assign(&coset_omega);
                    }
                });
            }
        });

        let mut final_values = vec![F::zero(); new_size];

        let results_ref = &results;

        worker.scope(final_values.len(), |scope, chunk| {
            for (i, v) in final_values.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = i * chunk;
                    for v in v.iter_mut() {
                        let coset_idx = idx % factor;
                        let element_idx = idx / factor;
                        *v = results_ref[coset_idx][element_idx];

                        idx += 1;
                    }
                });
            }
        });

        Polynomial::from_values(final_values)
    }

    // pub fn lde_using_multiple_cosets_with_precomputation<P: FftPrecomputations<F>>(
    //     self,
    //     worker: &Worker,
    //     factor: usize,
    //     precomputed_omegas: &P,
    // ) -> Result<Polynomial<F, Values>, SynthesisError> {
    //     debug_assert!(self.coeffs.len().is_power_of_two());
    //     debug_assert_eq!(self.size(), precomputed_omegas.domain_size());

    //     if factor == 1 {
    //         return Ok(self.fft(&worker));
    //     }

    //     let num_cpus = worker.cpus;
    //     let num_cpus_hint = if num_cpus <= factor {
    //         Some(1)
    //     } else {
    //         let threads_per_coset = factor / num_cpus;
    //         // let mut threads_per_coset = factor / num_cpus;
    //         // if factor % num_cpus != 0 {
    //         //     if (threads_per_coset + 1).is_power_of_two() {
    //         //         threads_per_coset += 1;
    //         //     }
    //         // }
    //         Some(threads_per_coset)
    //     };

    //     assert!(factor.is_power_of_two());
    //     let new_size = self.coeffs.len() * factor;
    //     let domain = Domain::<F>::new_for_size(new_size as u64)?;

    //     let mut results = vec![self.coeffs; factor];

    //     let coset_omega = domain.generator;
    //     let this_domain_omega = self.omega;

    //     let log_n = self.exp;

    //     worker.scope(results.len(), |scope, chunk| {
    //         for (i, r) in results.chunks_mut(chunk).enumerate() {
    //             scope.spawn(move |_| {
    //                 let mut coset_generator = coset_omega.pow(&[i as u64]);
    //                 for r in r.iter_mut() {
    //                     distribute_powers(&mut r[..], &worker, coset_generator);
    //                     with_precomputation::fft::best_fft(
    //                         &mut r[..],
    //                         &worker,
    //                         &this_domain_omega,
    //                         log_n,
    //                         num_cpus_hint,
    //                         precomputed_omegas,
    //                     );
    //                     coset_generator.mul_assign(&coset_omega);
    //                 }
    //             });
    //         }
    //     });

    //     let mut final_values = vec![F::zero(); new_size];

    //     let results_ref = &results;

    //     worker.scope(final_values.len(), |scope, chunk| {
    //         for (i, v) in final_values.chunks_mut(chunk).enumerate() {
    //             scope.spawn(move |_| {
    //                 let mut idx = i * chunk;
    //                 for v in v.iter_mut() {
    //                     let coset_idx = idx % factor;
    //                     let element_idx = idx / factor;
    //                     *v = results_ref[coset_idx][element_idx];

    //                     idx += 1;
    //                 }
    //             });
    //         }
    //     });

    //     Polynomial::from_values(final_values)
    // }

    pub fn lde_using_bitreversed_ntt<P: CTPrecomputations<F>>(
        self,
        worker: &Worker,
        factor: usize,
        precomputed_omegas: &P,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        use std::time::Instant;
        debug_assert!(self.coeffs.len().is_power_of_two());
        debug_assert_eq!(self.size(), precomputed_omegas.domain_size());

        if factor == 1 {
            return Ok(self.fft(&worker));
        }

        let num_cpus = worker.cpus;
        let num_cpus_hint = if num_cpus <= factor {
            Some(1)
        } else {
            let threads_per_coset = factor / num_cpus;
            // let mut threads_per_coset = factor / num_cpus;
            // if factor % num_cpus != 0 {
            //     if (threads_per_coset + 1).is_power_of_two() {
            //         threads_per_coset += 1;
            //     }
            // }
            Some(threads_per_coset)
        };

        assert!(factor.is_power_of_two());
        let current_size = self.coeffs.len();
        let new_size = current_size * factor;
        let domain = Domain::<F>::new_for_size(new_size as u64)?;

        let mut results = vec![self.coeffs; factor];

        let coset_omega = domain.generator;

        let log_n_u32 = self.exp;
        let log_n = log_n_u32 as usize;

        // for r in results.iter_mut().skip(1) {
        //     let mut coset_generator = coset_omega;
        //     distribute_powers(&mut r[..], &worker, coset_generator);
        //     coset_generator.mul_assign(&coset_omega);
        // }

        worker.scope(results.len(), |scope, chunk| {
            for (i, r) in results.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut coset_generator = coset_omega.pow(&[i as u64]);
                    let mut gen_power = i;
                    for r in r.iter_mut() {
                        if gen_power != 0 {
                            distribute_powers_serial(&mut r[..], coset_generator);
                        }
                        // distribute_powers(&mut r[..], &worker, coset_generator);
                        cooley_tukey_ntt::best_ct_ntt(
                            &mut r[..],
                            &worker,
                            log_n_u32,
                            num_cpus_hint,
                            precomputed_omegas,
                        );

                        coset_generator.mul_assign(&coset_omega);
                        gen_power += 1;
                    }
                });
            }
        });

        // let mut final_values = vec![F::zero(); new_size];

        let mut final_values = Vec::with_capacity(new_size);
        unsafe { final_values.set_len(new_size) };

        // copy here is more complicated: to have the value in a natural order
        // one has to use coset_idx to address the result element
        // and use bit-reversed lookup for an element index

        let results_ref = &results;

        worker.scope(final_values.len(), |scope, chunk| {
            for (i, v) in final_values.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = i * chunk;
                    for v in v.iter_mut() {
                        let coset_idx = idx % factor;
                        let element_idx = idx / factor;
                        let element_idx = cooley_tukey_ntt::bitreverse(element_idx, log_n);
                        *v = results_ref[coset_idx][element_idx];

                        idx += 1;
                    }
                });
            }
        });

        // let res_ptr = &mut final_values[..] as *mut [F];

        // let factor_log_n = crate::plonk::commitments::transparent::utils::log2_floor(factor);

        //  worker.scope(results.len(), |scope, chunk| {
        //     for (chunk_idx, r) in results.chunks(chunk).enumerate() {
        //         let res = unsafe {&mut *res_ptr};
        //         scope.spawn(move |_| {
        //             // elements from the coset i should be on the places
        //             // of sequence i, i + lde_factor, i + 2*lde_factor, ...
        //             let mut coset_idx = chunk_idx * chunk;
        //             for r in r.iter() {
        //                 for (element_idx, el) in r.iter().enumerate() {
        //                     let write_to = (cooley_tukey_ntt::bitreverse(element_idx, log_n) << factor_log_n) | coset_idx;
        //                     res[write_to] = *el;
        //                 }

        //                 coset_idx += 1;
        //             }
        //         });
        //     }
        // });

        Polynomial::from_values(final_values)
    }

    // pub fn lde_using_bitreversed_ntt_no_allocations_lowest_bits_reversed<P: CTPrecomputations<F>>(
    //     self,
    //     worker: &Worker,
    //     factor: usize,
    //     precomputed_omegas: &P
    // ) -> Result<Polynomial<F, Values>, SynthesisError> {
    //     debug_assert!(self.coeffs.len().is_power_of_two());
    //     debug_assert_eq!(self.size(), precomputed_omegas.domain_size());

    //     if factor == 1 {
    //         return Ok(self.fft(&worker));
    //     }

    //     let num_cpus = worker.cpus;
    //     let num_cpus_hint = if num_cpus <= factor {
    //         Some(1)
    //     } else {
    //         let threads_per_coset = factor / num_cpus;
    //         // let mut threads_per_coset = factor / num_cpus;
    //         // if factor % num_cpus != 0 {
    //         //     if (threads_per_coset + 1).is_power_of_two() {
    //         //         threads_per_coset += 1;
    //         //     }
    //         // }
    //         Some(threads_per_coset)
    //     };

    //     assert!(factor.is_power_of_two());
    //     let current_size = self.coeffs.len();
    //     let new_size = self.coeffs.len() * factor;
    //     let domain = Domain::<F>::new_for_size(new_size as u64)?;

    //     // let mut results = vec![self.coeffs.clone(); factor];

    //     let mut result = Vec::with_capacity(new_size);
    //     unsafe { result.set_len(new_size)};

    //     let r = &mut result[..] as *mut [F];

    //     let coset_omega = domain.generator;

    //     let log_n = self.exp;

    //     let range: Vec<usize> = (0..factor).collect();

    //     let self_coeffs_ref = &self.coeffs;

    //     // copy

    //     worker.scope(range.len(), |scope, chunk| {
    //         for coset_idx in range.chunks(chunk) {
    //             let r = unsafe {&mut *r};
    //             scope.spawn(move |_| {
    //                 for coset_idx in coset_idx.iter() {
    //                     let start = current_size * coset_idx;
    //                     let end = start + current_size;
    //                     let copy_start_pointer: *mut F = r[start..end].as_mut_ptr();

    //                     unsafe { std::ptr::copy_nonoverlapping(self_coeffs_ref.as_ptr(), copy_start_pointer, current_size) };
    //                 }
    //             });
    //         }
    //     });

    //     // let mut coset_generator = F::one();
    //     // for _ in 0..factor {
    //     //     result.extend_from_slice(&self.coeffs);
    //     //     // if coset_idx != 0 {
    //     //     //     let start = coset_idx * current_size;
    //     //     //     let end = start + current_size;
    //     //     //     distribute_powers(&mut result[start..end], &worker, coset_generator);
    //     //     //     // cooley_tukey_ntt::best_ct_ntt(&mut result[start..end], &worker, log_n, Some(log_num_cpus as usize), precomputed_omegas);
    //     //     // }
    //     //     // coset_generator.mul_assign(&coset_omega);
    //     // }
    //     // println!("Copying taken {:?}", start.elapsed());

    //     // for coset_idx in 0..factor {
    //     //     result.extend_from_slice(&self.coeffs);
    //     //     if coset_idx != 0 {
    //     //         let start = coset_idx * current_size;
    //     //         let end = start + current_size;
    //     //         distribute_powers(&mut result[start..end], &worker, coset_generator);
    //     //     }
    //     //     coset_generator.mul_assign(&coset_omega);
    //     // }

    //     // for r in results.iter_mut().skip(1) {
    //     //     let mut coset_generator = coset_omega;
    //     //     distribute_powers(&mut r[..], &worker, coset_generator);
    //     //     coset_generator.mul_assign(&coset_omega);
    //     // }

    //     // let start = Instant::now();

    //     let to_spawn = worker.cpus;

    //     let chunk = Worker::chunk_size_for_num_spawned_threads(factor, to_spawn);

    //     worker.scope(0, |scope, _| {
    //         for thread_id in 0..to_spawn {
    //             let r = unsafe {&mut *r};
    //             scope.spawn(move |_| {
    //                 let start = thread_id * chunk;
    //                 let end = if start + chunk <= factor {
    //                     start + chunk
    //                 } else {
    //                     factor
    //                 };
    //                 let mut gen_power = start;
    //                 let mut coset_generator = coset_omega.pow(&[start as u64]);
    //                 for i in start..end {
    //                     let from = current_size * i;
    //                     let to = from + current_size;
    //                     if gen_power != 0 {
    //                         distribute_powers_serial(&mut r[from..to], coset_generator);
    //                     }
    //                     cooley_tukey_ntt::best_ct_ntt(&mut r[from..to], &worker, log_n, num_cpus_hint, precomputed_omegas);
    //                     coset_generator.mul_assign(&coset_omega);
    //                     gen_power += 1;
    //                 }
    //             });
    //         }
    //     });

    //     // println!("NTT taken {:?}", start.elapsed());

    //     // let start = Instant::now();

    //     // let mut final_values = vec![F::zero(); new_size];

    //     // println!("Allocation of result taken {:?}", start.elapsed());

    //     // let results_ref = &results;

    //     // copy here is more complicated: to have the value in a natural order
    //     // one has to use coset_idx to address the result element
    //     // and use bit-reversed lookup for an element index

    //     // let log_n = log_n as usize;

    //     // let start = Instant::now();

    //     // let total_len = result.len();

    //     // let chunk = Worker::chunk_size_for_num_spawned_threads(total_len, to_spawn);

    //     // let lower_bits_mask = (1 << log_n) - 1;

    //     // let higher_bits_mask = !lower_bits_mask;

    //     // worker.scope(0, |scope, _| {
    //     //     for thread_id in 0..to_spawn {
    //     //         let r = unsafe {&mut *r};
    //     //         scope.spawn(move |_| {
    //     //             let start = thread_id * chunk;
    //     //             let end = if start + chunk <= total_len {
    //     //                 start + chunk
    //     //             } else {
    //     //                 total_len
    //     //             };
    //     //             for j in start..end {
    //     //                 let element_idx = j & lower_bits_mask;
    //     //                 let coset_mask = j & higher_bits_mask;
    //     //                 let rj = cooley_tukey_ntt::bitreverse(element_idx, log_n) | coset_mask;
    //     //                 if j < rj {
    //     //                     r.swap(j, rj);
    //     //                 }
    //     //             }
    //     //         });
    //     //     }
    //     // });

    //     // println!("Final copying taken {:?}", start.elapsed());

    //     Polynomial::from_values(result)
    // }

    pub fn coset_filtering_lde(
        mut self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(self.coset_fft(&worker));
        }
        assert!(factor.is_power_of_two());
        self.distribute_powers(worker, F::multiplicative_generator());

        let new_size = self.coeffs.len() * factor;
        let domain = Domain::new_for_size(new_size as u64)?;

        let mut lde = self.coeffs;
        lde.resize(new_size as usize, F::zero());
        best_lde(
            &mut lde,
            worker,
            &domain.generator,
            domain.power_of_two as u32,
            factor,
        );

        Polynomial::from_values(lde)
    }

    pub fn coset_lde_using_multiple_cosets_naive(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(self.coset_fft(&worker));
        }
        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::new_for_size(new_size as u64)?;

        let mut results = vec![];

        let mut coset_generator = F::multiplicative_generator();

        for _ in 0..factor {
            let coeffs = self.clone();
            let lde = coeffs.coset_fft_for_generator(&worker, coset_generator);

            results.push(lde.into_coeffs());
            coset_generator.mul_assign(&domain.generator);
        }

        let mut final_values = vec![F::zero(); new_size];

        let results_ref = &results;

        worker.scope(final_values.len(), |scope, chunk| {
            for (i, v) in final_values.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = i * chunk;
                    for v in v.iter_mut() {
                        let coset_idx = idx % factor;
                        let element_idx = idx / factor;
                        *v = results_ref[coset_idx][element_idx];

                        idx += 1;
                    }
                });
            }
        });

        Polynomial::from_values(final_values)
    }

    pub fn coset_lde_using_multiple_cosets(
        self,
        worker: &Worker,
        factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());

        if factor == 1 {
            return Ok(self.coset_fft(&worker));
        }

        let num_cpus = worker.cpus;
        let num_cpus_hint = if num_cpus <= factor {
            Some(1)
        } else {
            let threads_per_coset = factor / num_cpus;
            // let mut threads_per_coset = factor / num_cpus;
            // if factor % num_cpus != 0 {
            //     if (threads_per_coset + 1).is_power_of_two() {
            //         threads_per_coset += 1;
            //     }
            // }
            Some(threads_per_coset)
        };

        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::<F>::new_for_size(new_size as u64)?;

        let mut results = vec![self.coeffs; factor];

        let coset_omega = domain.generator;
        let this_domain_omega = self.omega;

        let log_n = self.exp;

        worker.scope(results.len(), |scope, chunk| {
            for (i, r) in results.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut coset_generator = coset_omega.pow(&[i as u64]);
                    coset_generator.mul_assign(&F::multiplicative_generator());
                    for r in r.iter_mut() {
                        distribute_powers_serial(&mut r[..], coset_generator);
                        // distribute_powers(&mut r[..], &worker, coset_generator);
                        best_fft(
                            &mut r[..],
                            &worker,
                            &this_domain_omega,
                            log_n,
                            // num_cpus_hint,
                            None,
                        );
                        coset_generator.mul_assign(&coset_omega);
                    }
                });
            }
        });

        let mut final_values = vec![F::zero(); new_size];

        let results_ref = &results;

        worker.scope(final_values.len(), |scope, chunk| {
            for (i, v) in final_values.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = i * chunk;
                    for v in v.iter_mut() {
                        let coset_idx = idx % factor;
                        let element_idx = idx / factor;
                        *v = results_ref[coset_idx][element_idx];

                        idx += 1;
                    }
                });
            }
        });

        Polynomial::from_values(final_values)
    }

    pub fn fft(mut self, worker: &Worker) -> Polynomial<F, Values> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        best_fft(&mut self.coeffs, worker, &self.omega, self.exp, None);

        Polynomial::<F, Values> {
            coeffs: self.coeffs,
            exp: self.exp,
            omega: self.omega,
            omegainv: self.omegainv,
            geninv: self.geninv,
            minv: self.minv,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn coset_fft(mut self, worker: &Worker) -> Polynomial<F, Values> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        self.distribute_powers(worker, F::multiplicative_generator());
        self.fft(worker)
    }

    pub fn coset_fft_for_generator(mut self, worker: &Worker, gen: F) -> Polynomial<F, Values> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        self.distribute_powers(worker, gen);
        self.fft(worker)
    }

    pub fn add_assign(&mut self, worker: &Worker, other: &Polynomial<F, Coefficients>) {
        assert!(self.coeffs.len() >= other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.add_assign(&b);
                    }
                });
            }
        });
    }

    pub fn add_assign_scaled(
        &mut self,
        worker: &Worker,
        other: &Polynomial<F, Coefficients>,
        scaling: &F,
    ) {
        assert!(self.coeffs.len() >= other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        let mut tmp = *b;
                        tmp.mul_assign(&scaling);
                        a.add_assign(&tmp);
                    }
                });
            }
        });
    }

    /// Perform O(n) subtraction of one polynomial from another in the domain.
    pub fn sub_assign(&mut self, worker: &Worker, other: &Polynomial<F, Coefficients>) {
        assert!(self.coeffs.len() >= other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.sub_assign(&b);
                    }
                });
            }
        });
    }

    pub fn sub_assign_scaled(
        &mut self,
        worker: &Worker,
        other: &Polynomial<F, Coefficients>,
        scaling: &F,
    ) {
        assert!(self.coeffs.len() >= other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        let mut tmp = *b;
                        tmp.mul_assign(&scaling);
                        a.sub_assign(&tmp);
                    }
                });
            }
        });
    }

    pub fn evaluate_at(&self, worker: &Worker, g: F) -> F {
        let num_threads = worker.get_num_spawned_threads(self.coeffs.len());
        let mut subvalues = vec![F::zero(); num_threads];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (i, (a, s)) in self
                .coeffs
                .chunks(chunk)
                .zip(subvalues.chunks_mut(1))
                .enumerate()
            {
                scope.spawn(move |_| {
                    let mut x = g.pow([(i * chunk) as u64]);
                    for a in a.iter() {
                        let mut value = x;
                        value.mul_assign(&a);
                        s[0].add_assign(&value);
                        x.mul_assign(&g);
                    }
                });
            }
        });

        let mut result = F::zero();
        for v in subvalues.iter() {
            result.add_assign(&v);
        }

        result
    }
}

impl<F: PrimeField> Polynomial<F, Values> {
    pub fn new_for_size(size: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        let coeffs = vec![F::zero(); size];

        Self::from_values(coeffs)
    }

    pub fn from_values(mut values: Vec<F>) -> Result<Polynomial<F, Values>, SynthesisError> {
        let coeffs_len = values.len();

        let domain = Domain::new_for_size(coeffs_len as u64)?;
        let exp = domain.power_of_two as u32;
        let m = domain.size as usize;
        let omega = domain.generator;

        values.resize(m, F::zero());

        Ok(Polynomial::<F, Values> {
            coeffs: values,
            exp: exp,
            omega: omega,
            omegainv: omega.inverse().unwrap(),
            geninv: F::multiplicative_generator().inverse().unwrap(),
            minv: F::from_str(&format!("{}", m)).unwrap().inverse().unwrap(),
            _marker: std::marker::PhantomData,
        })
    }

    pub fn from_values_unpadded(values: Vec<F>) -> Result<Polynomial<F, Values>, SynthesisError> {
        let coeffs_len = values.len();

        let domain = Domain::new_for_size(coeffs_len as u64)?;
        let exp = domain.power_of_two as u32;
        let m = domain.size as usize;
        let omega = domain.generator;

        Ok(Polynomial::<F, Values> {
            coeffs: values,
            exp: exp,
            omega: omega,
            omegainv: omega.inverse().unwrap(),
            geninv: F::multiplicative_generator().inverse().unwrap(),
            minv: F::from_str(&format!("{}", m)).unwrap().inverse().unwrap(),
            _marker: std::marker::PhantomData,
        })
    }

    // this function should only be used on the values that are bitreverse enumerated
    pub fn clone_subset_assuming_bitreversed(
        &self,
        subset_factor: usize,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        if subset_factor == 1 {
            return Ok(self.clone());
        }

        assert!(subset_factor.is_power_of_two());

        let current_size = self.coeffs.len();
        let new_size = current_size / subset_factor;

        let mut result = Vec::with_capacity(new_size);
        unsafe { result.set_len(new_size) };

        // copy elements. If factor is 2 then non-reversed we would output only elements that are == 0 mod 2
        // If factor is 2 and we are bit-reversed - we need to only output first half of the coefficients
        // If factor is 4 then we need to output only the first 4th part
        // if factor is 8 - only the first 8th part

        let start = 0;
        let end = new_size;

        result.copy_from_slice(&self.coeffs[start..end]);

        // unsafe { result.set_len(new_size)};
        // let copy_to_start_pointer: *mut F = result[..].as_mut_ptr();
        // let copy_from_start_pointer: *const F = self.coeffs[start..end].as_ptr();

        // unsafe { std::ptr::copy_nonoverlapping(copy_from_start_pointer, copy_to_start_pointer, new_size) };

        Polynomial::from_values(result)
    }

    pub fn pow(&mut self, worker: &Worker, exp: u64) {
        if exp == 2 {
            return self.square(&worker);
        }
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        *v = v.pow([exp]);
                    }
                });
            }
        });
    }

    pub fn square(&mut self, worker: &Worker) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        v.square();
                    }
                });
            }
        });
    }

    pub fn ifft(mut self, worker: &Worker) -> Polynomial<F, Coefficients> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        best_fft(&mut self.coeffs, worker, &self.omegainv, self.exp, None);

        worker.scope(self.coeffs.len(), |scope, chunk| {
            let minv = self.minv;

            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v {
                        v.mul_assign(&minv);
                    }
                });
            }
        });

        Polynomial::<F, Coefficients> {
            coeffs: self.coeffs,
            exp: self.exp,
            omega: self.omega,
            omegainv: self.omegainv,
            geninv: self.geninv,
            minv: self.minv,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn icoset_fft(self, worker: &Worker) -> Polynomial<F, Coefficients> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        let geninv = self.geninv;
        let mut res = self.ifft(worker);
        res.distribute_powers(worker, geninv);

        res
    }

    pub fn icoset_fft_for_generator(
        self,
        worker: &Worker,
        coset_generator: &F,
    ) -> Polynomial<F, Coefficients> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        let geninv = coset_generator.inverse().expect("must exist");
        let mut res = self.ifft(worker);
        res.distribute_powers(worker, geninv);

        res
    }

    pub fn add_assign(&mut self, worker: &Worker, other: &Polynomial<F, Values>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.add_assign(&b);
                    }
                });
            }
        });
    }

    pub fn add_constant(&mut self, worker: &Worker, constant: &F) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for a in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for a in a.iter_mut() {
                        a.add_assign(&constant);
                    }
                });
            }
        });
    }

    pub fn sub_constant(&mut self, worker: &Worker, constant: &F) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for a in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for a in a.iter_mut() {
                        a.sub_assign(&constant);
                    }
                });
            }
        });
    }

    pub fn rotate(mut self, by: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        let mut values: Vec<_> = self.coeffs.drain(by..).collect();

        for c in self.coeffs.into_iter() {
            values.push(c);
        }

        Polynomial::from_values(values)
    }

    pub fn barycentric_evaluate_at(&self, worker: &Worker, g: F) -> Result<F, SynthesisError> {
        // use a barycentric formula

        // L_i(X) = (omega^i / N) / (X - omega^i) * (X^N - 1)
        // we'll have to pay more for batch inversion at some point, but
        // it's still useful
        let domain_size = self.size() as u64;
        assert!(domain_size.is_power_of_two());

        let mut vanishing_at_g = g.pow(&[domain_size]);
        vanishing_at_g.sub_assign(&F::one());

        // now generate (X - omega^i)
        let mut tmp = vec![g; domain_size as usize];

        let generator = self.omega;

        // constant factor = 1 / ( (1 / N) * (X^N - 1) ) = N / (X^N - 1)

        worker.scope(tmp.len(), |scope, chunk| {
            for (i, vals) in tmp.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    // let mut one_over_omega_pow = generator_inv.pow([(i*chunk) as u64]);
                    // one_over_omega_pow.mul_assign(&constant_factor);
                    let mut omega_power = generator.pow([(i * chunk) as u64]);
                    for val in vals.iter_mut() {
                        val.sub_assign(&omega_power); // (X - omega^i)
                                                      // val.mul_assign(&one_over_omega_pow); // (X - omega^i) * N / (X^N - 1) * omega^(-i), so when we inverse it's valid evaluation
                        omega_power.mul_assign(&generator);
                        // one_over_omega_pow.mul_assign(&generator_inv);
                    }
                });
            }
        });

        let mut values = Polynomial::from_values(tmp)?;
        values.batch_inversion(&worker)?;

        // now multiply by omega^i / N * (X^N - 1) and value for L_i(X)

        let mut constant_factor = vanishing_at_g;
        constant_factor.mul_assign(&self.minv);

        worker.scope(values.size(), |scope, chunk| {
            for (i, (vals, coeffs)) in values
                .as_mut()
                .chunks_mut(chunk)
                .zip(self.coeffs.chunks(chunk))
                .enumerate()
            {
                scope.spawn(move |_| {
                    let mut omega_power = generator.pow([(i * chunk) as u64]);
                    omega_power.mul_assign(&constant_factor);
                    for (val, coeff) in vals.iter_mut().zip(coeffs.iter()) {
                        val.mul_assign(&omega_power);
                        val.mul_assign(coeff);
                        omega_power.mul_assign(&generator);
                    }
                });
            }
        });

        values.calculate_sum(&worker)
    }


    pub fn calculate_shifted_grand_product(
        &self,
        worker: &Worker,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        let mut result = vec![F::zero(); self.coeffs.len() + 1];
        result[0] = F::one();

        // let not_shifted_product = self.calculate_grand_product(&worker)?;
        // result[1..].copy_from_slice(&not_shifted_product.into_coeffs()[..]);

        // Polynomial::from_values_unpadded(result)

        let work_chunk = &mut result[1..];
        assert!(work_chunk.len() == self.coeffs.len());

        let num_threads = worker.get_num_spawned_threads(self.coeffs.len());
        let mut subproducts = vec![F::one(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for ((g, c), s) in work_chunk
                .chunks_mut(chunk)
                .zip(self.coeffs.chunks(chunk))
                .zip(subproducts.chunks_mut(1))
            {
                scope.spawn(move |_| {
                    for (g, c) in g.iter_mut().zip(c.iter()) {
                        s[0].mul_assign(&c);
                        *g = s[0];
                    }
                });
            }
        });

        // subproducts are [abc, def, xyz]

        // we do not need the last one
        subproducts.pop().expect("has at least one value");

        let mut tmp = F::one();
        for s in subproducts.iter_mut() {
            tmp.mul_assign(&s);
            *s = tmp;
        }

        let chunk_len = worker.get_chunk_size(self.coeffs.len());

        worker.scope(0, |scope, _| {
            for (g, s) in work_chunk[chunk_len..]
                .chunks_mut(chunk_len)
                .zip(subproducts.chunks(1))
            {
                scope.spawn(move |_| {
                    for g in g.iter_mut() {
                        g.mul_assign(&s[0]);
                    }
                });
            }
        });

        Polynomial::from_values(result)
    }

    pub fn calculate_grand_product(
        &self,
        worker: &Worker,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        let mut result = vec![F::zero(); self.coeffs.len()];

        let num_threads = worker.get_num_spawned_threads(self.coeffs.len());
        let mut subproducts = vec![F::one(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for ((g, c), s) in result
                .chunks_mut(chunk)
                .zip(self.coeffs.chunks(chunk))
                .zip(subproducts.chunks_mut(1))
            {
                scope.spawn(move |_| {
                    for (g, c) in g.iter_mut().zip(c.iter()) {
                        s[0].mul_assign(&c);
                        *g = s[0];
                    }
                });
            }
        });

        // subproducts are [abc, def, xyz]

        // we do not need the last one
        subproducts.pop().expect("has at least one value");

        let mut tmp = F::one();
        for s in subproducts.iter_mut() {
            tmp.mul_assign(&s);
            *s = tmp;
        }

        let chunk_len = worker.get_chunk_size(self.coeffs.len());

        worker.scope(0, |scope, _| {
            for (g, s) in result[chunk_len..]
                .chunks_mut(chunk_len)
                .zip(subproducts.chunks(1))
            {
                scope.spawn(move |_| {
                    let c = s[0];
                    for g in g.iter_mut() {
                        g.mul_assign(&c);
                    }
                });
            }
        });

        Polynomial::from_values_unpadded(result)
    }

    pub fn calculate_grand_product_serial(&self) -> Result<Polynomial<F, Values>, SynthesisError> {
        let mut result = Vec::with_capacity(self.coeffs.len());

        let mut tmp = F::one();
        for c in self.coeffs.iter() {
            tmp.mul_assign(&c);
            result.push(tmp);
        }

        Polynomial::from_values_unpadded(result)
    }

    pub fn calculate_sum(&self, worker: &Worker) -> Result<F, SynthesisError> {
        let num_threads = worker.get_num_spawned_threads(self.coeffs.len());
        let mut subresults = vec![F::zero(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (c, s) in self.coeffs.chunks(chunk).zip(subresults.chunks_mut(1)) {
                scope.spawn(move |_| {
                    for c in c.iter() {
                        s[0].add_assign(&c);
                    }
                });
            }
        });

        let mut sum = F::zero();

        for el in subresults.iter() {
            sum.add_assign(&el);
        }

        Ok(sum)
    }

    pub fn calculate_grand_sum(
        &self,
        worker: &Worker,
    ) -> Result<(F, Polynomial<F, Values>), SynthesisError> {
        // first value is zero, then first element, then first + second, ...
        let mut result = vec![F::zero(); self.coeffs.len() + 1];

        let num_threads = worker.get_num_spawned_threads(self.coeffs.len());
        let mut subsums = vec![F::zero(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for ((g, c), s) in result[1..]
                .chunks_mut(chunk)
                .zip(self.coeffs.chunks(chunk))
                .zip(subsums.chunks_mut(1))
            {
                scope.spawn(move |_| {
                    for (g, c) in g.iter_mut().zip(c.iter()) {
                        s[0].add_assign(&c);
                        *g = s[0];
                    }
                });
            }
        });

        // subsums are [a+b+c, d+e+f, x+y+z]

        let mut tmp = F::zero();
        for s in subsums.iter_mut() {
            tmp.add_assign(&s);
            *s = tmp;
        }

        // sum over the full domain is the last element
        let domain_sum = subsums.pop().expect("has at least one value");

        let chunk_len = worker.get_chunk_size(self.coeffs.len());

        worker.scope(0, |scope, _| {
            for (g, s) in result[(chunk_len + 1)..]
                .chunks_mut(chunk_len)
                .zip(subsums.chunks(1))
            {
                scope.spawn(move |_| {
                    let c = s[0];
                    for g in g.iter_mut() {
                        g.add_assign(&c);
                    }
                });
            }
        });

        // let result = result.drain(1..).collect();

        let alt_total_sum = result.pop().expect("must pop the last element");

        assert_eq!(alt_total_sum, domain_sum);

        Ok((domain_sum, Polynomial::from_values_unpadded(result)?))
    }

    pub fn add_assign_scaled(
        &mut self,
        worker: &Worker,
        other: &Polynomial<F, Values>,
        scaling: &F,
    ) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        let mut tmp = *b;
                        tmp.mul_assign(&scaling);
                        a.add_assign(&tmp);
                    }
                });
            }
        });
    }

    /// Perform O(n) subtraction of one polynomial from another in the domain.
    pub fn sub_assign(&mut self, worker: &Worker, other: &Polynomial<F, Values>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.sub_assign(&b);
                    }
                });
            }
        });
    }

    pub fn sub_assign_scaled(
        &mut self,
        worker: &Worker,
        other: &Polynomial<F, Values>,
        scaling: &F,
    ) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        let mut tmp = *b;
                        tmp.mul_assign(&scaling);
                        a.sub_assign(&tmp);
                    }
                });
            }
        });
    }

    /// Perform O(n) multiplication of two polynomials in the domain.
    pub fn mul_assign(&mut self, worker: &Worker, other: &Polynomial<F, Values>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.mul_assign(&b);
                    }
                });
            }
        });
    }

    pub fn batch_inversion(&mut self, worker: &Worker) -> Result<(), SynthesisError> {
        let num_threads = worker.get_num_spawned_threads(self.coeffs.len());

        let mut grand_products = vec![F::one(); self.coeffs.len()];
        let mut subproducts = vec![F::one(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for ((a, g), s) in self
                .coeffs
                .chunks(chunk)
                .zip(grand_products.chunks_mut(chunk))
                .zip(subproducts.chunks_mut(1))
            {
                scope.spawn(move |_| {
                    for (a, g) in a.iter().zip(g.iter_mut()) {
                        s[0].mul_assign(&a);
                        *g = s[0];
                    }
                });
            }
        });

        // now coeffs are [a, b, c, d, ..., z]
        // grand_products are [a, ab, abc, d, de, def, ...., xyz]
        // subproducts are [abc, def, xyz]
        // not guaranteed to have equal length

        let mut full_grand_product = F::one();
        for sub in subproducts.iter() {
            full_grand_product.mul_assign(sub);
        }

        let product_inverse = full_grand_product
            .inverse()
            .ok_or(SynthesisError::DivisionByZero("has no inverse".to_string()))?;

        // now let's get [abc^-1, def^-1, ..., xyz^-1];
        let mut subinverses = vec![F::one(); num_threads];
        for (i, s) in subinverses.iter_mut().enumerate() {
            let mut tmp = product_inverse;
            for (j, p) in subproducts.iter().enumerate() {
                if i == j {
                    continue;
                }
                tmp.mul_assign(&p);
            }

            *s = tmp;
        }

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for ((a, g), s) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(grand_products.chunks(chunk))
                .zip(subinverses.chunks_mut(1))
            {
                scope.spawn(move |_| {
                    for (a, g) in a
                        .iter_mut()
                        .rev()
                        .zip(g.iter().rev().skip(1).chain(Some(F::one()).iter()))
                    {
                        // s[0] = abc^-1
                        // a = c
                        // g = ab
                        let tmp = *a; // c
                        *a = *g;
                        a.mul_assign(&s[0]); // a = ab*(abc^-1) = c^-1
                        s[0].mul_assign(&tmp); // s[0] = (ab)^-1
                    }
                });
            }
        });

        Ok(())
    }

    pub fn pop_last(&mut self) -> Result<F, SynthesisError> {
        let last = self.coeffs.pop().ok_or(SynthesisError::AssignmentMissing)?;

        Ok(last)
    }

    pub fn clone_shifted_assuming_natural_ordering(
        &self,
        by: usize,
    ) -> Result<Self, SynthesisError> {
        let len = self.coeffs.len();
        assert!(by < len);
        let mut new_coeffs = vec![F::zero(); self.coeffs.len()];
        new_coeffs[..(len - by)].copy_from_slice(&self.coeffs[by..]);
        new_coeffs[(len - by)..].copy_from_slice(&self.coeffs[..by]);

        Self::from_values(new_coeffs)
    }

    pub fn clone_shifted_assuming_bitreversed(
        &self,
        by: usize,
        worker: &Worker,
    ) -> Result<Self, SynthesisError> {
        let len = self.coeffs.len();
        assert!(by < len);
        let mut extended_clone = Vec::with_capacity(len + by);
        extended_clone.extend_from_slice(&self.coeffs);
        let mut tmp = Self::from_values(extended_clone)?;
        tmp.bitreverse_enumeration(&worker);

        let mut coeffs = tmp.into_coeffs();
        let tmp: Vec<_> = coeffs.drain(..by).collect();
        coeffs.extend(tmp);

        let mut tmp = Self::from_values(coeffs)?;
        tmp.bitreverse_enumeration(&worker);

        Ok(tmp)
    }
}

// impl<F: PartialTwoBitReductionField> Polynomial<F, Coefficients> {
//     pub fn bitreversed_lde_using_bitreversed_ntt_with_partial_reduction<P: CTPrecomputations<F>>(
//         self,
//         worker: &Worker,
//         factor: usize,
//         precomputed_omegas: &P,
//         coset_factor: &F,
//     ) -> Result<Polynomial<F, Values>, SynthesisError> {
//         debug_assert!(self.coeffs.len().is_power_of_two());
//         debug_assert_eq!(self.size(), precomputed_omegas.domain_size());

//         if factor == 1 {
//             return Ok(self.fft(&worker));
//         }

//         let num_cpus = worker.cpus;
//         let num_cpus_hint = if num_cpus <= factor {
//             Some(1)
//         } else {
//             let mut threads_per_coset = num_cpus / factor;
//             if threads_per_coset == 0 {
//                 threads_per_coset = 1;
//             } else if num_cpus % factor != 0 {
//                 threads_per_coset += 1;
//             }
//             // let mut threads_per_coset = factor / num_cpus;
//             // if factor % num_cpus != 0 {
//             //     if (threads_per_coset + 1).is_power_of_two() {
//             //         threads_per_coset += 1;
//             //     }
//             // }
//             Some(threads_per_coset)
//         };

//         assert!(factor.is_power_of_two());
//         let current_size = self.coeffs.len();
//         let new_size = self.coeffs.len() * factor;
//         let domain = Domain::<F>::new_for_size(new_size as u64)?;

//         // let mut results = vec![self.coeffs.clone(); factor];

//         let mut result = Vec::with_capacity(new_size);
//         unsafe { result.set_len(new_size) };

//         let r = &mut result[..] as *mut [F];

//         let coset_omega = domain.generator;

//         let log_n = self.exp;

//         let range: Vec<usize> = (0..factor).collect();

//         let self_coeffs_ref = &self.coeffs;

//         // copy

//         worker.scope(range.len(), |scope, chunk| {
//             for coset_idx in range.chunks(chunk) {
//                 let r = unsafe { &mut *r };
//                 scope.spawn(move |_| {
//                     for coset_idx in coset_idx.iter() {
//                         let start = current_size * coset_idx;
//                         let end = start + current_size;
//                         let copy_start_pointer: *mut F = r[start..end].as_mut_ptr();

//                         unsafe {
//                             std::ptr::copy_nonoverlapping(
//                                 self_coeffs_ref.as_ptr(),
//                                 copy_start_pointer,
//                                 current_size,
//                             )
//                         };
//                     }
//                 });
//             }
//         });

//         let to_spawn = factor;
//         let coset_size = current_size;

//         use crate::plonk::commitments::transparent::utils::log2_floor;

//         let factor_log = log2_floor(factor) as usize;

//         // let chunk = Worker::chunk_size_for_num_spawned_threads(factor, to_spawn);

//         // Each coset will produce values at specific indexes only, e.g
//         // coset factor of omega^0 = 1 will produce elements that are only at places == 0 mod 16
//         // coset factor of omega^1 will produce elements that are only at places == 1 mod 16
//         // etc. We expect the output to be bitreversed, so
//         // elements for coset factor of omega^0 = 1 will need to be placed first (00 top bits, bitreversed 00)
//         // elements for coset factor of omega^1 will need to be placed after the first half (10 top bits, bitreversed 01)

//         worker.scope(0, |scope, _| {
//             for coset_idx in 0..to_spawn {
//                 let r = unsafe { &mut *r };
//                 scope.spawn(move |_| {
//                     let from = coset_size * coset_idx;
//                     let to = from + coset_size;
//                     let one = F::one();
//                     let bitreversed_power = cooley_tukey_ntt::bitreverse(coset_idx, factor_log);
//                     let mut coset_generator = coset_omega.pow(&[bitreversed_power as u64]);
//                     coset_generator.mul_assign(&coset_factor);
//                     if coset_generator != one {
//                         distribute_powers_with_num_cpus(
//                             &mut r[from..to],
//                             &worker,
//                             coset_generator,
//                             num_cpus_hint.expect("is some"),
//                         );
//                     }
//                     partial_reduction::best_ct_ntt_partial_reduction(
//                         &mut r[from..to],
//                         &worker,
//                         log_n,
//                         num_cpus_hint,
//                         precomputed_omegas,
//                     );
//                 });
//             }
//         });

//         Polynomial::from_values(result)
//     }
// }

// impl<F: PartialTwoBitReductionField> Polynomial<F, Values> {
//     pub fn ifft_using_bitreversed_ntt_with_partial_reduction<P: CTPrecomputations<F>>(
//         self,
//         worker: &Worker,
//         precomputed_omegas: &P,
//         coset_generator: &F,
//     ) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
//         if self.coeffs.len() <= worker.cpus * 4 {
//             return Ok(self.ifft(&worker));
//         }

//         let mut coeffs: Vec<_> = self.coeffs;
//         let exp = self.exp;
//         cooley_tukey_ntt::partial_reduction::best_ct_ntt_partial_reduction(
//             &mut coeffs,
//             worker,
//             exp,
//             Some(worker.cpus),
//             precomputed_omegas,
//         );

//         let mut this = Polynomial::from_coeffs(coeffs)?;

//         this.bitreverse_enumeration(&worker);

//         let geninv = coset_generator.inverse().expect("must exist");

//         worker.scope(this.coeffs.len(), |scope, chunk| {
//             let minv = this.minv;

//             for v in this.coeffs.chunks_mut(chunk) {
//                 scope.spawn(move |_| {
//                     for v in v {
//                         v.mul_assign(&minv);
//                     }
//                 });
//             }
//         });

//         if geninv != F::one() {
//             this.distribute_powers(&worker, geninv);
//         }

//         Ok(this)
//     }
// }

impl<F: PrimeField> Polynomial<F, Values> {
    /// taken in natural enumeration
    /// outputs in natural enumeration
    pub fn ifft_using_bitreversed_ntt<P: CTPrecomputations<F>>(
        self,
        worker: &Worker,
        precomputed_omegas: &P,
        coset_generator: &F,
    ) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
        if self.coeffs.len() <= worker.cpus * 4 {
            return Ok(self.ifft(&worker));
        }

        let mut coeffs: Vec<_> = self.coeffs;
        let exp = self.exp;
        cooley_tukey_ntt::best_ct_ntt(
            &mut coeffs,
            worker,
            exp,
            Some(worker.cpus),
            precomputed_omegas,
        );
        let mut this = Polynomial::from_coeffs(coeffs)?;

        this.bitreverse_enumeration(&worker);

        let geninv = coset_generator.inverse().expect("must exist");

        worker.scope(this.coeffs.len(), |scope, chunk| {
            let minv = this.minv;

            for v in this.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v {
                        v.mul_assign(&minv);
                    }
                });
            }
        });

        if geninv != F::one() {
            this.distribute_powers(&worker, geninv);
        }

        Ok(this)
    }
}

impl<F: PrimeField> Polynomial<F, Coefficients> {
    pub fn bitreversed_lde_using_bitreversed_ntt<P: CTPrecomputations<F>>(
        self,
        worker: &Worker,
        factor: usize,
        precomputed_omegas: &P,
        coset_factor: &F,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        debug_assert!(self.coeffs.len().is_power_of_two());
        debug_assert_eq!(self.size(), precomputed_omegas.domain_size());

        if factor == 1 {
            return Ok(self.fft(&worker));
        }

        let num_cpus = worker.cpus;
        let num_cpus_hint = if num_cpus <= factor {
            Some(1)
        } else {
            let mut threads_per_coset = num_cpus / factor;
            if threads_per_coset == 0 {
                threads_per_coset = 1;
            } else if num_cpus % factor != 0 {
                threads_per_coset += 1;
            }
            // let mut threads_per_coset = factor / num_cpus;
            // if factor % num_cpus != 0 {
            //     if (threads_per_coset + 1).is_power_of_two() {
            //         threads_per_coset += 1;
            //     }
            // }
            Some(threads_per_coset)
        };

        assert!(factor.is_power_of_two());
        let current_size = self.coeffs.len();
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::<F>::new_for_size(new_size as u64)?;

        // let mut results = vec![self.coeffs.clone(); factor];

        let mut result = Vec::with_capacity(new_size);
        unsafe { result.set_len(new_size) };

        let r = &mut result[..] as *mut [F];

        let coset_omega = domain.generator;

        let log_n = self.exp;

        let range: Vec<usize> = (0..factor).collect();

        let self_coeffs_ref = &self.coeffs;

        // copy

        worker.scope(range.len(), |scope, chunk| {
            for coset_idx in range.chunks(chunk) {
                let r = unsafe { &mut *r };
                scope.spawn(move |_| {
                    for coset_idx in coset_idx.iter() {
                        let start = current_size * coset_idx;
                        let end = start + current_size;
                        let copy_start_pointer: *mut F = r[start..end].as_mut_ptr();

                        unsafe {
                            std::ptr::copy_nonoverlapping(
                                self_coeffs_ref.as_ptr(),
                                copy_start_pointer,
                                current_size,
                            )
                        };
                    }
                });
            }
        });

        let to_spawn = factor;
        let coset_size = current_size;

        // use crate::plonk::commitments::transparent::utils::log2_floor;
        use crate::utils::log2_floor;

        let factor_log = log2_floor(factor) as usize;

        // let chunk = Worker::chunk_size_for_num_spawned_threads(factor, to_spawn);

        // Each coset will produce values at specific indexes only, e.g
        // coset factor of omega^0 = 1 will produce elements that are only at places == 0 mod 16
        // coset factor of omega^1 will produce elements that are only at places == 1 mod 16
        // etc. We expect the output to be bitreversed, so
        // elements for coset factor of omega^0 = 1 will need to be placed first (00 top bits, bitreversed 00)
        // elements for coset factor of omega^1 will need to be placed after the first half (10 top bits, bitreversed 01)

        worker.scope(0, |scope, _| {
            for coset_idx in 0..to_spawn {
                let r = unsafe { &mut *r };
                scope.spawn(move |_| {
                    let from = coset_size * coset_idx;
                    let to = from + coset_size;
                    let one = F::one();
                    let bitreversed_power = cooley_tukey_ntt::bitreverse(coset_idx, factor_log);
                    let mut coset_generator = coset_omega.pow(&[bitreversed_power as u64]);
                    coset_generator.mul_assign(&coset_factor);
                    if coset_generator != one {
                        distribute_powers_with_num_cpus(
                            &mut r[from..to],
                            &worker,
                            coset_generator,
                            num_cpus_hint.expect("is some"),
                        );
                    }
                    cooley_tukey_ntt::best_ct_ntt(
                        &mut r[from..to],
                        &worker,
                        log_n,
                        num_cpus_hint,
                        precomputed_omegas,
                    );
                });
            }
        });

        Polynomial::from_values(result)
    }

    /// taken in natural enumeration
    /// outputs in natural enumeration
    pub fn fft_using_bitreversed_ntt<P: CTPrecomputations<F>>(
        self,
        worker: &Worker,
        precomputed_omegas: &P,
        coset_generator: &F,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        if self.coeffs.len() <= worker.cpus * 4 {
            return Ok(self.coset_fft_for_generator(&worker, *coset_generator));
        }

        let mut this = self;
        if coset_generator != &F::one() {
            this.distribute_powers(&worker, *coset_generator);
        }

        let mut coeffs: Vec<_> = this.coeffs;
        let exp = this.exp;
        cooley_tukey_ntt::best_ct_ntt(
            &mut coeffs,
            worker,
            exp,
            Some(worker.cpus),
            precomputed_omegas,
        );
        let mut this = Polynomial::from_values(coeffs)?;

        this.bitreverse_enumeration(&worker);

        Ok(this)
    }

    /// taken in natural enumeration
    /// outputs in natural enumeration
    pub fn fft_using_bitreversed_ntt_output_bitreversed<P: CTPrecomputations<F>>(
        self,
        worker: &Worker,
        precomputed_omegas: &P,
        coset_generator: &F,
    ) -> Result<Polynomial<F, Values>, SynthesisError> {
        if self.coeffs.len() <= worker.cpus * 4 {
            return Ok(self.coset_fft_for_generator(&worker, *coset_generator));
        }

        let mut this = self;
        if coset_generator != &F::one() {
            this.distribute_powers(&worker, *coset_generator);
        }

        let mut coeffs: Vec<_> = this.coeffs;
        let exp = this.exp;
        cooley_tukey_ntt::best_ct_ntt(
            &mut coeffs,
            worker,
            exp,
            Some(worker.cpus),
            precomputed_omegas,
        );
        let this = Polynomial::from_values(coeffs)?;

        Ok(this)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::ff::{Field, PrimeField};
    use crate::Fr;
    use crate::fft::multicore::*;
    use crate::fft::cooley_tukey_ntt::*;
    use rand::Rng;

    fn make_random_field_elements<F: PrimeField>(
        worker: &Worker,
        num_elements: usize,
    ) -> Vec<F> {
        use rand::{XorShiftRng, SeedableRng, Rand, Rng, ChaChaRng};
    
        let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        make_random_field_elements_for_rng(worker, num_elements, rng)
    }

    fn make_random_field_elements_for_rng<F: PrimeField, R: Rng>(
        worker: &Worker,
        num_elements: usize,
        mut rng: R
    ) -> Vec<F> {
        let mut result = vec![F::zero(); num_elements];

        use rand::{XorShiftRng, SeedableRng, Rand, Rng, ChaChaRng};

        worker.scope(result.len(), |scope, chunk| {
            for r in result.chunks_mut(chunk)
            {
                let seed: [u32; 4] = rng.gen();
                let subrng = ChaChaRng::from_seed(&seed);
                scope.spawn(move |_| {
                    let mut subrng = subrng;
                    for r in r.iter_mut() {
                        *r = subrng.gen();
                    }
                });
            }
        });

        result 
    }

    #[test]
    fn compare_ntt_and_naive_lde() {
        let size: usize = 1 << 4;

        println!("size {} next power of two {}", size, size.next_power_of_two());

        let lde_factor = 2;

        let size_with_lde = size * lde_factor;

        let worker = Worker::new();
        

        let omegas_bitreversed =
            BitReversedOmegas::<Fr>::new_for_domain_size(size.next_power_of_two());

        let omegas_inv_bitreversed = <OmegasInvBitreversed::<Fr> as CTPrecomputations::<Fr>>::new_for_domain_size(size.next_power_of_two());

        let coeffs = make_random_field_elements::<Fr>(&worker, size);

        let poly = Polynomial::from_coeffs(coeffs).unwrap();

        let lde_original = poly.clone().lde(&worker, lde_factor).expect("");

        let lde_ntt = poly
            .clone()
            .bitreversed_lde_using_bitreversed_ntt(
                &worker,
                lde_factor,
                &omegas_bitreversed,
                &Fr::multiplicative_generator(),
            )
            .expect("");

        let mut expected_evals = lde_original.clone();
        expected_evals.mul_assign(&worker, &lde_original);

        let mut actual_evals = lde_ntt.clone();
        actual_evals.mul_assign(&worker, &lde_ntt);

        let expected = expected_evals.ifft(&worker);

        // let omegas_bitreversed_2n =
        //     BitReversedOmegas::<Fr>::new_for_domain_size(size_with_lde.next_power_of_two());

        let mut actual = actual_evals.ifft_using_bitreversed_ntt(
            &worker,
            // &omegas_bitreversed_2n,
            &omegas_inv_bitreversed,
            &Fr::multiplicative_generator(),
        ).expect("");

        actual.bitreverse_enumeration(&worker);

        assert_eq!(expected, actual);
    }
    #[test]
    fn test_addition_comparison_between_ntt_and_naive_fft() {

        let size: usize = 1 << 4;

        println!("size {} next power of two {}", size, size.next_power_of_two());

        let lde_factor = 2;

        let size_with_lde = size * lde_factor;

        let worker = Worker::new();


        let omegas_bitreversed =
            BitReversedOmegas::<Fr>::new_for_domain_size(size.next_power_of_two());

        let omegas_inv_bitreversed = <OmegasInvBitreversed::<Fr> as CTPrecomputations::<Fr>>::new_for_domain_size(size.next_power_of_two());

        let coeffs = make_random_field_elements::<Fr>(&worker, size);

        let poly = Polynomial::from_coeffs(coeffs).unwrap();

        let original_evals = poly.clone().fft(&worker);

        let mut expected_evals = original_evals.clone();
        expected_evals.add_assign(&worker, &original_evals);

        let expected = expected_evals.ifft(&worker);
        
        // let coset_factor = &Fr::multiplicative_generator();
        let coset_factor = &Fr::one();

        // outputs natural enumeration so we don't need to change permutation
        let ntt_evals = poly
            .clone()
            .fft_using_bitreversed_ntt(
                &worker,
                &omegas_bitreversed,
                &coset_factor,
            )
            .expect("");

        let mut actual_evals = ntt_evals.clone();
        actual_evals.add_assign(&worker, &ntt_evals);

        let mut actual = actual_evals.icoset_fft_for_generator(
            &worker,
            &coset_factor,
        );

        assert_eq!(expected.as_ref(), actual.as_ref());
    }
    #[test]
    fn test_multiplication_comparison_between_ntt_and_naive_lde() {
        let size: usize = 1 << 4;

        println!("size {} next power of two {}", size, size.next_power_of_two());

        let lde_factor = 2;

        let size_with_lde = size * lde_factor;

        let worker = Worker::new();

        let omegas_bitreversed =
            BitReversedOmegas::<Fr>::new_for_domain_size(size.next_power_of_two());

        let omegas_inv_bitreversed = <OmegasInvBitreversed::<Fr> as CTPrecomputations::<Fr>>::new_for_domain_size(size.next_power_of_two());

        let coeffs = make_random_field_elements::<Fr>(&worker, size);

        let poly = Polynomial::from_coeffs(coeffs).unwrap();

        let original_evals = poly.clone().lde(&worker, lde_factor).expect("some lde");

        let mut expected_evals = original_evals.clone();
        expected_evals.mul_assign(&worker, &original_evals);

        let expected = expected_evals.ifft(&worker);
        
        let coset_factor = &Fr::multiplicative_generator();

        // outputs in bitreverse enumeration so we need to change permutation
        let ntt_evals = poly
            .clone()
            .bitreversed_lde_using_bitreversed_ntt(
                &worker,
                lde_factor,
                &omegas_bitreversed,
                &coset_factor,
            )
            .expect("some lde");

        let mut actual_evals = ntt_evals.clone();
        actual_evals.mul_assign(&worker, &ntt_evals);

        actual_evals.bitreverse_enumeration(&worker);
        let mut actual = actual_evals.icoset_fft_for_generator(
            &worker,
            &coset_factor,
        );

        assert_eq!(expected.as_ref(), actual.as_ref());
    }

    // #[test]
    // fn test_multiplication_comparison_between_ntt_and_naive_ifft() {        
    //     let size: usize = 1 << 4;

    //     println!("size {} next power of two {}", size, size.next_power_of_two());

    //     let lde_factor = 2;

    //     let size_with_lde = size * lde_factor;

    //     let worker = Worker::new();

    //     let omegas_bitreversed =
    //         BitReversedOmegas::<Fr>::new_for_domain_size(size.next_power_of_two());

    //     let omegas_inv_bitreversed = <OmegasInvBitreversed::<Fr> as CTPrecomputations::<Fr>>::new_for_domain_size(size.next_power_of_two());

    //     let coeffs = crate::kate_commitment::test::make_random_field_elements::<Fr>(&worker, size);

    //     let poly = Polynomial::from_coeffs(coeffs).unwrap();

    //     let original_evals = poly.clone().lde(&worker, lde_factor).expect("some lde");

    //     let mut expected_evals = original_evals.clone();
    //     expected_evals.mul_assign(&worker, &original_evals);

    //     let expected = expected_evals.ifft(&worker);
        
    //     let coset_factor = &Fr::multiplicative_generator();
    //     // let coset_factor = &Fr::one();

    //     // outputs in bitreverse enumeration so we need to change permutation
    //     let ntt_evals = poly
    //         .clone()
    //         .bitreversed_lde_using_bitreversed_ntt(
    //             &worker,
    //             lde_factor,
    //             &omegas_bitreversed,
    //             &coset_factor,
    //         )
    //         .expect("some lde");

    //     let mut actual_evals = ntt_evals.clone();
    //     actual_evals.mul_assign(&worker, &ntt_evals);

    //     actual_evals.bitreverse_enumeration(&worker);
    //     let mut actual = actual_evals.ifft_using_bitreversed_ntt(
    //         &worker,
    //         &omegas_inv_bitreversed,
    //         &coset_factor,
    //     ).expect("some poly");

    //     assert_eq!(expected.as_ref(), actual.as_ref());
    // }
}
