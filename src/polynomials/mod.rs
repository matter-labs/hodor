use ff::PrimeField;
use crate::domains::*;
use crate::SynthesisError;
use crate::fft::multicore::*;
use crate::fft::*;

pub trait PolynomialForm: Sized + Copy + Clone {}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Coefficients { }

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Values { }

impl PolynomialForm for Coefficients {}
impl PolynomialForm for Values{}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Polynomial<F: PrimeField, P: PolynomialForm> {
    coeffs: Vec<F>,
    pub exp: u32,
    pub omega: F,
    pub omegainv: F,
    pub geninv: F,
    pub minv: F,
    _marker: std::marker::PhantomData<P>
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

    pub fn distribute_powers(&mut self, worker: &Worker, g: F)
    {
        distribute_powers(&mut self.coeffs, &worker, g);
    }

    pub fn scale(&mut self, worker: &Worker, g: F)
    {
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

    pub fn negate(&mut self, worker: &Worker)
    {
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

    pub fn pad_by_factor(&mut self, factor: usize) -> Result<(), SynthesisError> {
        if factor == 1 {
            return Ok(());
        }
        let next_power_of_two = factor.next_power_of_two();
        if factor != next_power_of_two {
            return Err(SynthesisError::Error);
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
            return Err(SynthesisError::Error);
        }
        let next_power_of_two = new_size.next_power_of_two();
        if new_size != next_power_of_two {
            return Err(SynthesisError::Error);
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

    pub fn from_coeffs(mut coeffs: Vec<F>) -> Result<Polynomial<F, Coefficients>, SynthesisError>
    {
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
            _marker: std::marker::PhantomData
        })
    }

    pub fn from_roots(roots: Vec<F>, worker: &Worker) -> Result<Polynomial<F, Coefficients>, SynthesisError>
    {

        let coeffs_len = roots.len() + 1;

        let domain = Domain::<F>::new_for_size(coeffs_len as u64)?;
        let num_threads = worker.cpus;

        // vector of vectors of polynomial coefficients for subproducts
        let mut subterms = vec![vec![]; num_threads];

        worker.scope(roots.len(), |scope, chunk| {
            for (r, s) in roots.chunks(chunk)
                    .zip(subterms.chunks_mut(1)) {
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
        domain_size: u64
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
        domain_size: u64
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
    pub fn lde(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        self.lde_using_multiple_cosets(worker, factor)
        // self.filtering_lde(worker, factor)
    }

    #[inline(always)]
    pub fn coset_lde(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        self.coset_lde_using_multiple_cosets(worker, factor)
        // self.filtering_coset_lde(worker, factor)
    }


    pub fn filtering_lde(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        if factor == 1 {
            return Ok(self.fft(&worker));
        }
        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::new_for_size(new_size as u64)?;

        let mut lde = self.coeffs;
        lde.resize(new_size as usize, F::zero());
        best_lde(&mut lde, worker, &domain.generator, domain.power_of_two as u32, factor);

        Polynomial::from_values(lde)
    }

    pub fn lde_using_multiple_cosets_naive(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
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
                    let mut idx = i*chunk;
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

    pub fn lde_using_multiple_cosets(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        if factor == 1 {
            return Ok(self.fft(&worker));
        }

        let num_cpus = worker.cpus;
        let num_cpus_hint = if num_cpus <= factor {
            Some(1)
        } else {
            let mut threads_per_coset = factor / num_cpus;
            if factor % num_cpus != 0 {
                threads_per_coset += 1;
            }
            Some(threads_per_coset)
        };

        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::<F>::new_for_size(new_size as u64)?;

        let mut results = vec![vec![]; factor];

        let coset_omega = domain.generator;
        let this_domain_omega = self.omega;

        let coeffs = self.coeffs;
        let log_n = self.exp;

        worker.scope(results.len(), |scope, chunk| {
            for (i, r) in results.chunks_mut(chunk).enumerate() {
                let coeffs_for_coset = coeffs.clone();
                scope.spawn(move |_| {
                    let mut coset_generator = coset_omega.pow(&[i as u64]);
                    for r in r.iter_mut() {
                        let mut c = coeffs_for_coset.clone();
                        distribute_powers(&mut c, &worker, coset_generator);
                        best_fft(&mut c, &worker, &this_domain_omega, log_n, num_cpus_hint);
                        *r = c;
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
                    let mut idx = i*chunk;
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

    pub fn coset_filtering_lde(mut self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        if factor == 1 {
            return Ok(self.coset_fft(&worker));
        }
        assert!(factor.is_power_of_two());
        self.distribute_powers(worker, F::multiplicative_generator());

        let new_size = self.coeffs.len() * factor;
        let domain = Domain::new_for_size(new_size as u64)?;

        let mut lde = self.coeffs;
        lde.resize(new_size as usize, F::zero());
        best_lde(&mut lde, worker, &domain.generator, domain.power_of_two as u32, factor);

        Polynomial::from_values(lde)
    }

    pub fn coset_lde_using_multiple_cosets_naive(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
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
                    let mut idx = i*chunk;
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

    pub fn coset_lde_using_multiple_cosets(self, worker: &Worker, factor: usize) -> Result<Polynomial<F, Values>, SynthesisError> {
        if factor == 1 {
            return Ok(self.coset_fft(&worker));
        }

        let num_cpus = worker.cpus;
        let num_cpus_hint = if num_cpus <= factor {
            Some(1)
        } else {
            let mut threads_per_coset = factor / num_cpus;
            if factor % num_cpus != 0 {
                threads_per_coset += 1;
            }
            Some(threads_per_coset)
        };

        assert!(factor.is_power_of_two());
        let new_size = self.coeffs.len() * factor;
        let domain = Domain::<F>::new_for_size(new_size as u64)?;

        let mut results = vec![vec![]; factor];

        let coset_omega = domain.generator;
        let this_domain_omega = self.omega;

        let coeffs = self.coeffs;
        let log_n = self.exp;

        worker.scope(results.len(), |scope, chunk| {
            for (i, r) in results.chunks_mut(chunk).enumerate() {
                let coeffs_for_coset = coeffs.clone();
                scope.spawn(move |_| {
                    let mut coset_generator = coset_omega.pow(&[i as u64]);
                    coset_generator.mul_assign(&F::multiplicative_generator());
                    for r in r.iter_mut() {
                        let mut c = coeffs_for_coset.clone();
                        distribute_powers(&mut c, &worker, coset_generator);
                        best_fft(&mut c, &worker, &this_domain_omega, log_n, num_cpus_hint);
                        *r = c;
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
                    let mut idx = i*chunk;
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

    pub fn fft(mut self, worker: &Worker) -> Polynomial<F, Values>
    {
        best_fft(&mut self.coeffs, worker, &self.omega, self.exp, None);

        Polynomial::<F, Values> {
            coeffs: self.coeffs,
            exp: self.exp,
            omega: self.omega,
            omegainv: self.omegainv,
            geninv: self.geninv,
            minv: self.minv,
            _marker: std::marker::PhantomData
        }
    }

    pub fn coset_fft(mut self, worker: &Worker) -> Polynomial<F, Values>
    {
        self.distribute_powers(worker, F::multiplicative_generator());

        self.fft(worker)
    }

    pub fn coset_fft_for_generator(mut self, worker: &Worker, gen: F) -> Polynomial<F, Values>
    {
        self.distribute_powers(worker, gen);

        self.fft(worker)
    }

    pub fn add_assign(&mut self, worker: &Worker, other: &Polynomial<F, Coefficients>) {
        assert!(self.coeffs.len() >= other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.add_assign(&b);
                    }
                });
            }
        });
    }

    pub fn add_assign_scaled(&mut self, worker: &Worker, other: &Polynomial<F, Coefficients>, scaling: &F) {
        assert!(self.coeffs.len() >= other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
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
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.sub_assign(&b);
                    }
                });
            }
        });
    }

    pub fn evaluate_at(&self, worker: &Worker, g: F) -> F {
        let num_threads = worker.cpus;
        let mut subvalues = vec![F::zero(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (i, (a, s)) in self.coeffs.chunks(chunk)
                        .zip(subvalues.chunks_mut(1))
                        .enumerate() {
                scope.spawn(move |_| {
                    let mut x = g.pow([(i*chunk) as u64]);
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

    pub fn from_values(mut values: Vec<F>) -> Result<Polynomial<F, Values>, SynthesisError>
    {
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
            _marker: std::marker::PhantomData
        })
    }

    pub fn pow(&mut self, worker: &Worker, exp: u64)
    {
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

    pub fn square(&mut self, worker: &Worker)
    {
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

    pub fn ifft(mut self, worker: &Worker) -> Polynomial<F, Coefficients>
    {
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
            _marker: std::marker::PhantomData
        }
    }

    pub fn icoset_fft(self, worker: &Worker) -> Polynomial<F, Coefficients>
    {
        let geninv = self.geninv;
        let mut res = self.ifft(worker);
        res.distribute_powers(worker, geninv);

        res
    }

    pub fn icoset_fft_for_generator(self, worker: &Worker, geninv: F) -> Polynomial<F, Coefficients>
    {
        let mut res = self.ifft(worker);
        res.distribute_powers(worker, geninv);

        res
    }

    pub fn add_assign(&mut self, worker: &Worker, other: &Polynomial<F, Values>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
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

    pub fn add_assign_scaled(&mut self, worker: &Worker, other: &Polynomial<F, Values>, scaling: &F) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(other.coeffs.len(), |scope, chunk| {
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
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
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.sub_assign(&b);
                    }
                });
            }
        });
    }

    /// Perform O(n) multiplication of two polynomials in the domain.
    pub fn mul_assign(&mut self, worker: &Worker, other: &Polynomial<F, Values>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self.coeffs.chunks_mut(chunk).zip(other.coeffs.chunks(chunk)) {
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
            for ((a, g), s) in self.coeffs.chunks(chunk)
                        .zip(grand_products.chunks_mut(chunk))
                        .zip(subproducts.chunks_mut(1)) {
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

        let product_inverse = full_grand_product.inverse().ok_or(SynthesisError::Error)?;

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
            for ((a, g), s) in self.coeffs.chunks_mut(chunk)
                        .zip(grand_products.chunks(chunk))
                        .zip(subinverses.chunks_mut(1)) {
                scope.spawn(move |_| {
                    for (a, g) in a.iter_mut().rev()
                                .zip(g.iter().rev().skip(1).chain(Some(F::one()).iter())) {
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
}


#[test]
fn test_batch_inversion() {
    use crate::Fr;
    use crate::fft::multicore::*;
    use crate::ff::Field;

    for size in 1..256 {
        if !(size as usize).is_power_of_two() {
            continue;
        }
        println!("Size = {}", size);
        let mut inputs = vec![];
        for i in 1..=size {
            let f = Fr::from_str(&i.to_string()).unwrap();
            inputs.push(f);
        }

        let worker = Worker::new();

        let mut p = Polynomial::from_values(inputs.clone()).unwrap();
        p.batch_inversion(&worker).unwrap();

        for (k, (a, b)) in inputs.iter().zip(p.as_ref().iter()).enumerate() {
            let inv = a.inverse().unwrap();
            assert!(*b == inv, "invalid at {}", k);
        }
    }
}

#[test]
fn test_lde_correctness() {
    use rand::{XorShiftRng, SeedableRng, Rand};
    const LOG_N: usize = 2;
    const BASE: usize = 1 << LOG_N;
    const LOG_LDE: usize = 4;
    const LDE_FACTOR: usize = 1 << LOG_LDE;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    use ff::Field;
    use crate::experiments::Fr;
    use crate::fft::multicore::Worker;
    use crate::polynomials::Polynomial;
    use std::time::Instant;

    let worker = Worker::new();

    let mut coeffs = (0..BASE).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

    let poly = Polynomial::from_coeffs(coeffs.clone()).unwrap();

    let p0 = poly.clone();
    let now = Instant::now();
    let coset_lde = p0.lde_using_multiple_cosets(&worker, LDE_FACTOR).unwrap();
    println!("LDE with multiple cosets taken {}ms", now.elapsed().as_millis());

    let p1 = poly.clone();
    let now = Instant::now();
    let filtering_lde = p1.filtering_lde(&worker, LDE_FACTOR).unwrap();
    println!("filtering LDE taken {}ms", now.elapsed().as_millis());

    let c = coeffs.clone();
    coeffs.resize(BASE * LDE_FACTOR, Fr::zero());

    let p2 = Polynomial::from_coeffs(coeffs).unwrap();
    let now = Instant::now();
    let naive_lde = p2.fft(&worker);
    println!("Naive FFT taken {}ms", now.elapsed().as_millis());

    let f = filtering_lde.into_coeffs();
    let n = naive_lde.into_coeffs();
    let c = coset_lde.into_coeffs();

    assert!(f == n);
    assert!(c == n);
}


#[test]
fn test_coset_lde_correctness() {
    use rand::{XorShiftRng, SeedableRng, Rand};
    const LOG_N: usize = 2;
    const BASE: usize = 1 << LOG_N;
    const LOG_LDE: usize = 4;
    const LDE_FACTOR: usize = 1 << LOG_LDE;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    use ff::Field;
    use crate::experiments::Fr;
    use crate::fft::multicore::Worker;
    use crate::polynomials::Polynomial;
    use std::time::Instant;

    let worker = Worker::new();

    let mut coeffs = (0..BASE).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

    let poly = Polynomial::from_coeffs(coeffs.clone()).unwrap();

    let p0 = poly.clone();
    let now = Instant::now();
    // let coset_lde = p0.coset_lde(&worker, LDE_FACTOR).unwrap();
    let coset_lde = p0.coset_lde_using_multiple_cosets(&worker, LDE_FACTOR).unwrap();
    println!("LDE with multiple cosets taken {}ms", now.elapsed().as_millis());

    let p1 = poly.clone();
    let now = Instant::now();
    let filtering_lde = p1.coset_filtering_lde(&worker, LDE_FACTOR).unwrap();
    println!("filtering LDE taken {}ms", now.elapsed().as_millis());

    let c = coeffs.clone();
    coeffs.resize(BASE * LDE_FACTOR, Fr::zero());

    let p2 = Polynomial::from_coeffs(coeffs).unwrap();
    let now = Instant::now();
    let naive_lde = p2.coset_fft(&worker);
    println!("Naive FFT taken {}ms", now.elapsed().as_millis());

    let f = filtering_lde.into_coeffs();
    let n = naive_lde.into_coeffs();
    let c = coset_lde.into_coeffs();

    assert!(f == n);
    assert!(c == n);
}

#[test]
fn test_various_ldes() {
    use rand::{XorShiftRng, SeedableRng, Rand};
    const LOG_N: usize = 22;
    const BASE: usize = 1 << LOG_N;
    const LOG_LDE: usize = 4;
    const LDE_FACTOR: usize = 1 << LOG_LDE;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    use ff::Field;
    use crate::experiments::Fr;
    use crate::fft::multicore::Worker;
    use crate::polynomials::Polynomial;
    use std::time::Instant;

    let worker = Worker::new();

    let mut coeffs = (0..BASE).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

    let poly = Polynomial::from_coeffs(coeffs.clone()).unwrap();

    let p0 = poly.clone();
    let now = Instant::now();
    let coset_lde = p0.lde_using_multiple_cosets(&worker, LDE_FACTOR).unwrap();
    println!("LDE with multiple cosets taken {}ms", now.elapsed().as_millis());

    let p1 = poly.clone();
    let now = Instant::now();
    let filtering_lde = p1.filtering_lde(&worker, LDE_FACTOR).unwrap();
    println!("filtering LDE taken {}ms", now.elapsed().as_millis());

    let mut c = coeffs.clone();
    coeffs.resize(BASE * LDE_FACTOR, Fr::zero());

    let p2 = Polynomial::from_coeffs(coeffs).unwrap();
    let now = Instant::now();
    let naive_lde = p2.fft(&worker);
    println!("Naive FFT taken {}ms", now.elapsed().as_millis());

    let now = Instant::now();
    let _small_fft = poly.fft(&worker);
    println!("Small FFT taken {}ms", now.elapsed().as_millis());

    let now = Instant::now();
    crate::fft::fft::serial_fft(&mut c, &naive_lde.omega, LOG_N as u32);
    println!("Serial small FFT taken {}ms", now.elapsed().as_millis());

    let f = filtering_lde.into_coeffs();
    let n = naive_lde.into_coeffs();
    let c = coset_lde.into_coeffs();
    assert!(f == n);
    assert!(c == n);

}