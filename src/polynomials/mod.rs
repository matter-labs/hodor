use ff::PrimeField;
use crate::domains::*;
use crate::SynthesisError;
use crate::fft::multicore::*;
use crate::fft::fft::*;

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
    exp: u32,
    pub omega: F,
    omegainv: F,
    geninv: F,
    minv: F,
    _marker: std::marker::PhantomData<P>
}


impl<F: PrimeField, P: PolynomialForm> Polynomial<F, P> {
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
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (i, v) in self.coeffs.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = g.pow(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        v.mul_assign(&u);
                        u.mul_assign(&g);
                    }
                });
            }
        });
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

    // // TODO: implement fancier FFT and implement separate LDE functionality
    // pub fn extend(&mut self, factor: usize, worker: &Worker) -> Result<(), SynthesisError> {
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

    pub fn extend(&mut self, factor: usize, _worker: &Worker) -> Result<(), SynthesisError> {
        if factor == 1 {
            return Ok(());
        }
        let next_power_of_two = factor.next_power_of_two();
        if factor != next_power_of_two {
            return Err(SynthesisError::Error);
        }
        
        let new_size = self.coeffs.len() * factor;
        self.coeffs.resize(new_size, F::zero());

        Ok(())
    }

    pub fn fft(mut self, worker: &Worker) -> Polynomial<F, Values>
    {
        best_fft(&mut self.coeffs, worker, &self.omega, self.exp);

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

    /// Perform O(n) subtraction of one polynomial from another in the domain.
    pub fn sub_assign(&mut self, worker: &Worker, other: &Polynomial<F, Coefficients>) {
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

    pub fn evaluate_at(&mut self, worker: &Worker, g: F) -> F {
        let num_threads = worker.cpus;
        let mut subvalues = vec![F::zero(); num_threads as usize];

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (i, (a, s)) in self.coeffs.chunks(chunk)
                        .zip(subvalues.chunks_mut(chunk))
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

    pub fn ifft(mut self, worker: &Worker) -> Polynomial<F, Coefficients>
    {
        best_fft(&mut self.coeffs, worker, &self.omegainv, self.exp);

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

    pub fn batch_inversion(&mut self, worker: &Worker) {
        let num_threads = worker.cpus;
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

        let mut full_grand_product = F::one();
        for sub in subproducts.iter() {
            full_grand_product.mul_assign(sub);
        }

        let product_inverse = full_grand_product.inverse().expect("is non-zero");

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
    }
}
