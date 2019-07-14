use ff::PrimeField;
use crate::domains::*;
use crate::SynthesisError;
use crate::fft::multicore::*;
use crate::fft::fft::*;

pub trait PolynomialForm: Sized {}

#[derive(Debug)]
pub enum Coefficients { }

#[derive(Debug)]
pub enum Values { }

impl PolynomialForm for Coefficients {}
impl PolynomialForm for Values{}

#[derive(Debug)]
pub struct Polynomial<F: PrimeField, P: PolynomialForm> {
    coeffs: Vec<F>,
    exp: u32,
    omega: F,
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
}

impl<F: PrimeField> Polynomial<F, Coefficients> {
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

    /// This evaluates t(tau) for this domain, which is
    /// tau^m - 1 for these radix-2 domains.
    pub fn z(&self, tau: &F) -> F {
        let mut tmp = tau.pow(&[self.coeffs.len() as u64]);
        tmp.sub_assign(&F::one());

        tmp
    }

    /// The target polynomial is the zero polynomial in our
    /// evaluation domain, so we must perform division over
    /// a coset.
    pub fn divide_by_z_on_coset(&mut self, worker: &Worker)
    {
        let i = self.z(&F::multiplicative_generator()).inverse().unwrap();

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v {
                        v.mul_assign(&i);
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
}


impl<F: PrimeField> Polynomial<F, Values> {
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

}
