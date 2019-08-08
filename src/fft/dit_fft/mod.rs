use ff::PrimeField;
use super::multicore::*;

pub(crate) fn serial_dit_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32, non_zero_entries_count: usize)
{
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
    
    let n = a.len() as u64;
    assert_eq!(n, 1 << log_n);

    let mut m = 1;
    for _ in 0..log_n
    {
        let w_m = omega.pow(&[(m) as u64]);
        let block_len = n / m;

        for block in 0..m
        {
            let mut w = F::one();
            for k in (block * block_len)..(block * block_len + std::cmp::min(block_len/2, non_zero_entries_count as u64))
            {
                let t = a[(k + block_len / 2) as usize];
                let mut tmp = a[(k) as usize];
                tmp.sub_assign(&t);
                a[(k+ block_len / 2) as usize] = tmp;
                a[(k+ block_len / 2) as usize].mul_assign(&w);
                a[(k) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }
        }

        m *= 2;
    }

    for k in 0..n
    {
        let rk = bitreverse(k, log_n);
        if k < rk
        {
            a.swap(rk as usize, k as usize);
        }
    }
}

pub(crate) fn parallel_dit_fft<F: PrimeField>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    log_cpus: u32,
    non_zero_entries_count: usize
)
{
    assert!(log_n >= log_cpus);

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
                        let idx = i + (s << log_new_n);
                        let mut t = a[idx];
                        t.mul_assign(&elt);
                        tmp[i].add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_dit_fft(tmp, &new_omega, log_new_n, std::cmp::min(non_zero_entries_count, 1 << log_new_n));
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

pub(crate) fn best_dit_fft<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, non_zero_entries_count: usize)
{
    let log_cpus = worker.log_num_cpus();

    if log_n <= log_cpus {
        serial_dit_fft(a, omega, log_n, non_zero_entries_count);
    } else {
        parallel_dit_fft(a, worker, omega, log_n, log_cpus, non_zero_entries_count);
    }
}

