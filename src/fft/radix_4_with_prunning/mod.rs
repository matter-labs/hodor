use ff::PrimeField;
use super::multicore::*;


pub(crate) fn best_fft<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, non_zero_entries_count: usize)
{
    assert!(log_n % 2 == 0); 
    let mut log_cpus = worker.log_num_cpus();
    if (log_cpus % 2 != 0)
    {
        log_cpus -= 1;
    }

    // we split into radix-4 kernels, so we need more points to start
    if log_n <= log_cpus {
        serial_fft(a, omega, log_n, non_zero_entries_count);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus, non_zero_entries_count);
    }
}

#[inline(always)]
fn base_4_digit_reverse(mut n: u64, l: u64) -> u64 {
    let mut r = 0u64;
    for _ in 0..l {
        r = (r << 2) | (n & 3);
        n >>= 2;
    }

    r
}

#[inline(always)]
    fn bitreverse(mut n: u64, l: u32) -> u64
    {
        let mut r = 0;
        for _ in 0..l
        {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

pub(crate) fn serial_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32, non_zero_entries_count: usize)
{
    let n = a.len() as u64;
    assert_eq!(n, 1 << log_n);

    assert!(log_n % 2 == 0);

    // v = W_4
    let v = omega.pow(&[(n / 4) as u64]);

    let mut m = 1;
    for _ in 0..(log_n / 2)
    {
        let w_m = omega.pow(&[(m) as u64]);
        let mut block_len = n / (m);

        for block in 0..m
        {
            let mut w = F::one();
            for k in (block * block_len)..(block * block_len + std::cmp::min(block_len/4, (non_zero_entries_count as u64)))
            {
                // y_0 = x_0 + x_1 + x_2 + x_3
                // y_1 = x_0 + W_4 * x_1 - x_2 - W_4 * x_3
                // y_2 = x_0 - x_1 + x_2 - x3
                // y_3 = x_0 - W_4 * x_1 - x_2 + W_4 * x_3
                              
                let mut x0 = a[(k) as usize];
                let mut x1 = a[(k+block_len/4) as usize];
                let mut x2 = a[(k+block_len/2) as usize];
                let mut x3 = a[(k+block_len *3/4) as usize];
                
                let mut temp1 = x0;
                temp1.add_assign(&x2);

                let mut temp2 = x1;
                temp2.add_assign(&x3);

                a[(k) as usize] = temp1;
                a[(k) as usize].add_assign(&temp2);
                a[(k+block_len/2) as usize] = temp1;
                a[(k+block_len/2) as usize].sub_assign(&temp2);

                x1.mul_assign(&v);
                x3.mul_assign(&v);

                temp1 = x0;
                temp1.sub_assign(&x2);
                temp2 = x1;
                temp2.sub_assign(&x3);

                a[(k+block_len/4) as usize] = temp1;
                a[(k+block_len/4) as usize].add_assign(&temp2);
                a[(k+3*block_len/4) as usize] = temp1;
                a[(k+3*block_len/4) as usize].sub_assign(&temp2);

                let mut u = w;
                a[(k+block_len/4) as usize].mul_assign(&u);
                u.mul_assign(&w);
                a[(k+block_len/2) as usize].mul_assign(&u);
                u.mul_assign(&w);
                a[(k+3*block_len/4) as usize].mul_assign(&u);

                w.mul_assign(&w_m);
            }
        }

        m *= 4;
    }
    
    let num_digits = (log_n/2) as u64;
    for k in 0..n {
        let rk = base_4_digit_reverse(k, num_digits);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }
}

pub(crate) fn parallel_fft<F: PrimeField>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    log_cpus: u32,
    non_zero_entries_count: usize
)
{
    assert!(log_n >= log_cpus);
    
    //we need log_n and log_cpu to be even
    assert!(log_n % 2 == 0);
    assert!(log_cpus % 2 == 0);
    
    let num_cpus = 1 << log_cpus;
    let log_new_n = (log_n - log_cpus);
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
                        let idx = (i + (s << log_new_n));
                        let mut t = a[idx];
                        t.mul_assign(&elt);
                        tmp[i].add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n, std::cmp::min(non_zero_entries_count, 1 << log_new_n));
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
}
