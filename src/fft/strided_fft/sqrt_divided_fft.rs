
use ff::PrimeField;
use super::utils::*;
use super::super::multicore::*;
use super::super::fft::{subview_mut};
use super::shuffle::transpose;
use super::fft::small_size_serial_fft;

pub fn radix_sqrt<F: PrimeField, const OUTER_SIZE: usize, const INNER_SIZE: usize, const MAX_LOOP_UNROLL: usize>(
    values: &mut [F], 
    omega: &F,
    precomputed_twiddle_factors: &[F],
    worker: &Worker
)
{
    if values.len() <= 1 {
        return;
    }

    // Recurse by splitting along the square root
    // Round such that outer is larger.
    let length = values.len();
    assert!(length.is_power_of_two());

    assert_eq!(precomputed_twiddle_factors.len() * 2, OUTER_SIZE);

    assert!(OUTER_SIZE >= INNER_SIZE);
    let stretch = OUTER_SIZE / INNER_SIZE;

    debug_assert_eq!(omega.pow(&[values.len() as u64]), F::one());
    debug_assert!(OUTER_SIZE == INNER_SIZE || OUTER_SIZE == 2 * INNER_SIZE);
    debug_assert_eq!(OUTER_SIZE * INNER_SIZE, length);

    // shuffle 
    transpose(values, INNER_SIZE, stretch);

    {
        let values_view = subview_mut::<F, OUTER_SIZE>(values);
        debug_assert_eq!(values_view.len(), INNER_SIZE);
        worker.scope(values_view.len(),  |scope, chunk_size| {
            for subwork in values_view.chunks_mut(chunk_size) {
                scope.spawn(move |_| {
                    for s in subwork.iter_mut() {
                        let mut ss = *s;
                        small_size_serial_fft::<F, OUTER_SIZE, MAX_LOOP_UNROLL>(&mut ss, precomputed_twiddle_factors, 0, stretch, stretch);
                        *s = ss;
                    }
                });
            }
        });
    }

    // shuffle back
    transpose(values, INNER_SIZE, stretch);

    {
        let values_view = subview_mut::<F, OUTER_SIZE>(values);
        debug_assert_eq!(values_view.len(), INNER_SIZE);
        worker.scope(values_view.len(),  |scope, chunk_size| {
            for (chunk_idx, subwork) in values_view.chunks_mut(chunk_size).enumerate() {
                scope.spawn(move |_| {
                    let start = chunk_idx * chunk_size;
                    for (i, s) in subwork.iter_mut().enumerate() {
                        let mut ss = *s;
                        let idx = start + i;
                        if idx > 0 {
                            let i = bitreverse_index_for_constant_size::<INNER_SIZE>(i);
                            let inner_twiddle = omega.pow(&[i as u64]);
                            let mut outer_twiddle = inner_twiddle;
                            for element in ss.iter_mut().skip(1) {
                                element.mul_assign(&outer_twiddle);
                                outer_twiddle.mul_assign(&inner_twiddle);
                            }
                        }
                        small_size_serial_fft::<F, OUTER_SIZE, MAX_LOOP_UNROLL>(&mut ss, precomputed_twiddle_factors, 0, 1, 1);
                        *s = ss;
                    }
                });
            }
        });
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_sqrt_fft_strategy() {
        use crate::ff::{Field, PrimeField};
        use crate::experiments::fields::vdf_128_ef::Fr;
        use crate::domains::Domain;
        use crate::fft::multicore::Worker;
        use crate::fft::fft::*;
        use std::time::Instant;
        use super::super::utils::precompute_twiddle_factors_parallelized;

        use rand::{XorShiftRng, SeedableRng, Rand};
        let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        // this looks much faster for 2^20 case and slow for 2^24 case for now

        const LOG_N: usize = 24;
        const N: usize = 1 << LOG_N;

        const INNER_SIZE: usize = 1 << (LOG_N / 2);
        const OUTER_SIZE: usize = N / INNER_SIZE;

        let reference = (0..N).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    
        let worker = Worker::new();

        let domain = Domain::<Fr>::new_for_size(N as u64).unwrap();
        let omega = domain.generator;

        let log_cpus = worker.log_num_cpus();

        // let start = Instant::now();
        // crate::fft::fft::parallel_fft(&mut reference[..], &worker, &omega, LOG_N as u32, log_cpus);
        // println!("default parallel fft taken {:?}", start.elapsed());

        // const CACHE_32MB: usize = (1 << 25) / std::mem::size_of::<Fr>(); // 16MB cache
        // const CACHE_16MB: usize = (1 << 24) / std::mem::size_of::<Fr>(); // 16MB cache

        // let start = Instant::now();
        // parallel_cache_friendly_fft::<Fr, CACHE_32MB>(&mut input, &worker, &omega, LOG_N as u32);
        // println!("cache friendly fft taken {:?} for {} element cache", start.elapsed(), CACHE_32MB);

        // let start = Instant::now();
        // const K32: usize = N / CACHE_32MB;
        // parallel_partitioned_fft::<Fr, CACHE_32MB, K32>(&mut input, &worker, &omega, LOG_N as u32);
        // println!("hard partitioned fft taken {:?} for {} element cache", start.elapsed(), CACHE_16MB);

        // let start = Instant::now();
        // const K16: usize = N / CACHE_16MB;
        // parallel_partitioned_fft::<Fr, CACHE_16MB, K16>(&mut input, &worker, &omega, LOG_N as u32);
        // println!("hard partitioned fft taken {:?} for {} element cache", start.elapsed(), CACHE_16MB);

        // let start = Instant::now();
        // super::radix_sqrt::<_, OUTER_SIZE, INNER_SIZE, 128>(&mut input, &omega, &worker);
        // println!("SQRT strategy fft taken {:?} for outer = {}, inner = {}, loop unroll = {}", start.elapsed(), OUTER_SIZE, INNER_SIZE, 128);



        const INNER_SIZE_0: usize = 1 << (LOG_N / 2);
        const OUTER_SIZE_0: usize = N / INNER_SIZE_0;
        const UNROLL_0: usize = 128;

        let mut input_0 = reference.clone();

        let outer_omega = omega.pow(&[INNER_SIZE_0 as u64]);
        let twiddles = precompute_twiddle_factors_parallelized(&outer_omega, OUTER_SIZE, &worker);

        let start = Instant::now();
        super::radix_sqrt::<_, OUTER_SIZE_0, INNER_SIZE_0, UNROLL_0>(&mut input_0, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for outer = {}, inner = {}, loop unroll = {}", start.elapsed(), OUTER_SIZE_0, INNER_SIZE_0, UNROLL_0);

        // const INNER_SIZE_1: usize = 1 << 11;
        // const OUTER_SIZE_1: usize = N / INNER_SIZE_1;
        const UNROLL_1: usize = 256;

        let mut input = reference.clone();

        let start = Instant::now();
        super::radix_sqrt::<_, OUTER_SIZE_0, INNER_SIZE_0, UNROLL_1>(&mut input, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for outer = {}, inner = {}, loop unroll = {}", start.elapsed(), OUTER_SIZE_0, INNER_SIZE_0, UNROLL_1);

        assert_eq!(&input, &input_0);
        const UNROLL_2: usize = 512;

        let mut input = reference.clone();

        let start = Instant::now();
        super::radix_sqrt::<_, OUTER_SIZE_0, INNER_SIZE_0, UNROLL_2>(&mut input, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for outer = {}, inner = {}, loop unroll = {}", start.elapsed(), OUTER_SIZE_0, INNER_SIZE_0, UNROLL_2);

        assert_eq!(&input, &input_0);
        const UNROLL_3: usize = 1024;

        let mut input = reference.clone();

        let start = Instant::now();
        super::radix_sqrt::<_, OUTER_SIZE_0, INNER_SIZE_0, UNROLL_3>(&mut input, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for outer = {}, inner = {}, loop unroll = {}", start.elapsed(), OUTER_SIZE_0, INNER_SIZE_0, UNROLL_3);

        assert_eq!(&input, &input_0);


        // const INNER_SIZE_2: usize = 1 << 10;
        // const OUTER_SIZE_2: usize = N / INNER_SIZE_2;
        // const UNROLL_2: usize = 128;

        // let start = Instant::now();
        // super::radix_sqrt::<_, OUTER_SIZE_2, INNER_SIZE_2, UNROLL_2>(&mut input, &omega, &worker);
        // println!("SQRT strategy fft taken {:?} for outer = {}, inner = {}, loop unroll = {}", start.elapsed(), OUTER_SIZE_2, INNER_SIZE_2, UNROLL_2);
    }
}