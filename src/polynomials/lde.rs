use ff::PrimeField;
use super::multicore::*;

pub struct LDE<F: PrimeField> {
    coeffs: Vec<F>,
    lde_factor: usize,
    final_domain_len: usize,
    omega: F,
    zero: F
}

impl<F: PrimeField> LDE<F> {
    #[inline(always)]
    fn get(&self, idx: usize) -> &F {
        if idx < self.coeffs.len() {
            &self.coeffs[idx]
        } else if idx < self.final_domain_len {
            &self.zero
        } else {
            unreachable!();
        }
    }

    fn evaluate(self) -> Vec<F> {
        vec![]
    }

    pub(crate) fn best_eval(
        self, 
        worker: &Worker, 
    ) -> Vec<F>
    {
        let log_cpus = worker.log_num_cpus();
        let log_n = log2_floor(self.coeffs.len());

        if log_n <= log_cpus {
            Self::recursive_eval(&self.coeffs, self.omega, log_n)
        } else {
            self.parallel_eval(worker, log_n, log_cpus)
        }
    }


    pub(crate) fn recursive_eval(a: &[F], omega: F, log_n: u32) -> Vec<F>
    {
        let n = a.len();
        if n == 2 {
            let mut t0 = a[0];
            let mut t1 = t0;
            t0.add_assign(&a[1]);
            t1.sub_assign(&a[1]);

            return vec![t0, t1];
        }
        assert_eq!(n, 1 << log_n);

        let n_half = n / 2;

        // copy into subproblems
        let mut even = Vec::with_capacity(n_half);
        let mut odd = Vec::with_capacity(n_half);

        for c in a.chunks(2) {
            even.push(c[0]);
            odd.push(c[1]);
        }

        let mut w_n = omega;
        w_n.square();

        let next_log_n = log_n - 1;

        let y_0 = recursive_fft(&even, w_n, next_log_n);
        let y_1 = recursive_fft(&odd, w_n, next_log_n);

        let mut result = vec![F::zero(); n];

        let mut w = F::one();
        for (i, (y0, y1)) in y_0.into_iter()
                            .zip(y_1.into_iter())
                            .enumerate() {
            let mut tmp = y1;
            tmp.mul_assign(&w);

            result[i] = y0;
            result[i].add_assign(&tmp);

            result[i+n_half] = y0;
            result[i+n_half].sub_assign(&tmp);

            w.mul_assign(&omega);
        }

        result
    }


    pub(crate) fn parallel_eval(
        self,
        worker: &Worker,
        log_n: u32,
        log_cpus: u32
    ) -> Vec<F>
    {
        assert!(log_n >= log_cpus);

        let mut a = self.coeffs;

        let num_cpus = 1 << log_cpus;
        let log_new_n = log_n - log_cpus;
        let mut tmp = vec![vec![F::zero(); 1 << log_new_n]; num_cpus];
        let new_omega = omega.pow(&[num_cpus as u64]);

        worker.scope(0, |scope, _| {
            let a = &*a;

            for (j, tmp) in tmp.iter_mut().enumerate() {
                scope.spawn(move |_| {
                    // Shuffle into a sub-FFT
                    let omega_j = omega.pow(&[j as u64]);
                    let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                    let mut elt = F::one();
                    for i in 0..(1 << log_new_n) {
                        for s in 0..num_cpus {
                            let idx = (i + (s << log_new_n)) % (1 << log_n);
                            let mut t = a[idx];
                            t.mul_assign(&elt);
                            tmp[i].add_assign(&t);
                            elt.mul_assign(&omega_step);
                        }
                        elt.mul_assign(&omega_j);
                    }

                    // Perform sub-FFT
                    *tmp = recursive_fft(tmp, new_omega, log_new_n);
                });
            }
        });

        worker.scope(a.len(), |scope, chunk| {
            let tmp = &tmp;

            for (idx, a) in a.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = idx * chunk;
                    let mask = (1 << log_cpus) - 1;
                    for a in a {
                        *a = tmp[idx & mask][idx >> log_cpus];
                        idx += 1;
                    }
                });
            }
        });

        a
    }
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow+1)) <= num {
        pow += 1;
    }

    pow
}