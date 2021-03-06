use ff::PrimeField;
use super::utils::{bitreverse_index};
use super::super::multicore::*;
use super::shuffle::transpose;


#[inline(always)]
pub fn radix_2_butterfly<F: PrimeField>(values: &mut [F], offset: usize, stride: usize)
{
    // a + b, a - b
    unsafe {
        let i = offset;
        let j = offset + stride;
        let i_el = *values.get_unchecked(i);
        let j_el = *values.get_unchecked(j);
        values.get_unchecked_mut(i).add_assign(&j_el);
        *values.get_unchecked_mut(j) = i_el;
        values.get_unchecked_mut(j).sub_assign(&j_el);
    }
}

#[inline(always)]
pub fn radix_2_butterfly_with_twiddle<F: PrimeField>(values: &mut [F], twiddle: &F, offset: usize, stride: usize)
{
    // a + w*b, a - w*b

    // we can make use of partial reduction here:
    // a + w*b \in [0, 3p)
    // a + p - w*b \in [0, 2p)
    unsafe {
        let i = offset;
        let j = offset + stride;
        let i_el = *values.get_unchecked(i);
        let mut j_el = *values.get_unchecked(j);
        j_el.mul_assign(&twiddle);
        values.get_unchecked_mut(i).add_assign(&j_el);
        *values.get_unchecked_mut(j) = i_el;
        values.get_unchecked_mut(j).sub_assign(&j_el);
    }
}

pub fn non_generic_small_size_serial_fft<F: PrimeField, const MAX_LOOP_UNROLL: usize>(
    values: &mut [F],
    precomputed_twiddle_factors: &[F],
    offset: usize,
    count: usize,
    stride: usize,
) {
    debug_assert_eq!(values.len() % stride, 0);
    // work size
    let size = values.len() / stride;
    debug_assert!(size.is_power_of_two());
    debug_assert!(offset < stride);
    if size > 1 {
        // Inner FFT radix size/2 without explicit splitting
        if stride == count && count < MAX_LOOP_UNROLL {
            non_generic_small_size_serial_fft::<F, MAX_LOOP_UNROLL>(values, precomputed_twiddle_factors, offset, 2 * count, 2 * stride);
        } else {
            // we may parallelize this too as indexes do not overlap
            non_generic_small_size_serial_fft::<F, MAX_LOOP_UNROLL>(values, precomputed_twiddle_factors, offset, count, 2 * stride);
            non_generic_small_size_serial_fft::<F, MAX_LOOP_UNROLL>(values, precomputed_twiddle_factors, offset + stride, count, 2 * stride);
        }

        // unrolled loops
        // we expect them to be small enough in case of square root division strategy,
        // so we do not need to care about prefetches
        for i in offset..offset + count {
            radix_2_butterfly(values, i, stride);
        }
        for (offset, twiddle) in (offset..offset + size * stride)
            .step_by(2 * stride)
            .zip(precomputed_twiddle_factors.iter())
            .skip(1)
        {
            for i in offset..offset + count {
                radix_2_butterfly_with_twiddle(values, twiddle, i, stride)
            }
        }
    }
}

pub fn calcualate_inner_and_outer_sizes(size: usize) -> (usize, usize) {
    assert!(size.is_power_of_two());
    let log_n = size.trailing_zeros();
    let inner_size = 1 << (log_n / 2);
    let outer_size = size / inner_size; 


    (inner_size, outer_size)
}

pub fn non_generic_radix_sqrt<F: PrimeField, const MAX_LOOP_UNROLL: usize>(
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

    let (inner_size, outer_size) = calcualate_inner_and_outer_sizes(length);
    // TODO
    // assert_eq!(precomputed_twiddle_factors.len() * 2, outer_size);
    let stretch = outer_size / inner_size;

    debug_assert_eq!(omega.pow(&[values.len() as u64]), F::one());
    debug_assert!(outer_size == inner_size || outer_size == 2 * inner_size);
    debug_assert_eq!(outer_size * inner_size, length);

    // shuffle 
    transpose(values, inner_size, stretch);

    {
        worker.scope(inner_size,  |scope, num_inner_works_chunk| {
            // we parallelize inner_size units of work over M CPUs, each unit is of "outer size", so
            // split input values into chunks of outer_size * num_inner_works_chunk
            let mut full_slice = &mut *values;
            let mut num_spawned = (worker.num_cpus() as usize) / num_inner_works_chunk;
            if (worker.num_cpus() as usize) % inner_size != 0 {
                num_spawned += 1;
            }
            for _ in 0..num_spawned {
                let num_values = if num_inner_works_chunk*outer_size <= full_slice.len() {
                    num_inner_works_chunk*outer_size
                } else {
                    debug_assert_eq!(full_slice.len() % outer_size, 0);
                    full_slice.len()
                };
                let (subwork, rest) = full_slice.split_at_mut(num_values);
                full_slice = rest;
                scope.spawn(move |_| {
                    for s in subwork.chunks_mut(outer_size) {
                        non_generic_small_size_serial_fft::<F, MAX_LOOP_UNROLL>(s, precomputed_twiddle_factors, 0, stretch, stretch);
                    }
                });
            }
        });
    }

    // shuffle back
    transpose(values, inner_size, stretch);

    {
        worker.scope(inner_size,  |scope, num_inner_works_chunk| {
            let mut full_slice = values;
            let mut num_spawned = (worker.num_cpus() as usize) / num_inner_works_chunk;
            if (worker.num_cpus() as usize) % inner_size != 0 {
                num_spawned += 1;
            }
            for thread_idx in 0..num_spawned {
                let num_values = if num_inner_works_chunk*outer_size <= full_slice.len() {
                    num_inner_works_chunk*outer_size
                } else {
                    debug_assert_eq!(full_slice.len() % outer_size, 0);
                    full_slice.len()
                };
                let (subwork, rest) = full_slice.split_at_mut(num_values);
                full_slice = rest;
                scope.spawn(move |_| {
                    let start = thread_idx * num_inner_works_chunk;
                    for (i, s) in subwork.chunks_mut(outer_size).enumerate() {
                        let idx = start + i;
                        if idx > 0 {
                            let i = bitreverse_index(inner_size, i);
                            let inner_twiddle = omega.pow(&[i as u64]);
                            let mut outer_twiddle = inner_twiddle;
                            for element in s.iter_mut().skip(1) {
                                element.mul_assign(&outer_twiddle);
                                outer_twiddle.mul_assign(&inner_twiddle);
                            }
                        }
                        non_generic_small_size_serial_fft::<F, MAX_LOOP_UNROLL>(s, precomputed_twiddle_factors, 0, 1, 1);
                    }
                });
            }
        });
    }
}

#[cfg(test)]
mod test {

    use crate::{ff::{Field, PrimeField}, fft::fft::best_fft};
    use crate::experiments::fields::vdf_128_ef::Fr;
    use crate::domains::Domain;
    use crate::fft::multicore::Worker;
    use std::time::Instant;
    use super::{super::utils::precompute_twiddle_factors_parallelized, non_generic_radix_sqrt};
    use rand::{XorShiftRng, SeedableRng, Rand};

    #[test]
    fn test_nongeneric_sqrt_fft_strategy() {
 

        
        let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        const LOG_N: usize = 24;
        const N: usize = 1 << LOG_N;
    
        let inner_size = 1 << (LOG_N/2);
        let outer_size = 1 << (LOG_N - (LOG_N/2));

        let reference = (0..N).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    
        let worker = Worker::new();

        let domain = Domain::<Fr>::new_for_size(N as u64).unwrap();
        let omega = domain.generator;

        const UNROLL_0: usize = 128;

        let mut input_0 = reference.clone();

        let outer_omega = omega.pow(&[inner_size as u64]);
        let twiddles = precompute_twiddle_factors_parallelized(&outer_omega, outer_size, &worker);

        let start = Instant::now();
        super::non_generic_radix_sqrt::<_, UNROLL_0>(&mut input_0, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for loop unroll = {}", start.elapsed(), UNROLL_0);

        const UNROLL_1: usize = 256;

        let mut input = reference.clone();

        let start = Instant::now();
        super::non_generic_radix_sqrt::<_, UNROLL_1>(&mut input, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for loop unroll = {}", start.elapsed(), UNROLL_1);

        assert_eq!(&input, &input_0);
        const UNROLL_2: usize = 512;

        let mut input = reference.clone();

        let start = Instant::now();
        super::non_generic_radix_sqrt::<_, UNROLL_2>(&mut input, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for loop unroll = {}", start.elapsed(), UNROLL_2);

        assert_eq!(&input, &input_0);
        const UNROLL_3: usize = 1024;

        let mut input = reference.clone();

        let start = Instant::now();
        super::non_generic_radix_sqrt::<_, UNROLL_3>(&mut input, &omega, &twiddles, &worker);
        println!("SQRT strategy fft taken {:?} for loop unroll = {}", start.elapsed(), UNROLL_3);

        assert_eq!(&input, &input_0);
    }

    #[test]
    fn test_correctness_of_square_root_fft_comparison_against_best_fft(){
        let worker = Worker::new();
        let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let size: usize = 16; // 32*32

        let log_size  = size.trailing_zeros();

        assert_eq!(log_size, 4);

        let domain = Domain::new_for_size(size as u64).expect("a domain");

        let omega = domain.generator;

        let mut values: Vec<Fr> = std::iter::repeat(size).take(size).map(|_| Fr::rand(rng)).collect();

        let mut values_for_sqrt = values.clone();

        best_fft(&mut values, &worker, &omega, log_size, None);

        // square root fft is based on a vector of NxN values
        // so we should compute square root of size
        let inner_size = 1<<(log_size / 2);

        // we have a domain generator which w^(domain_size) = 1
        // so 
        // we need a inner_size root of unity  => w^(1<<inner_size) = 1
        let inner_generator = omega.pow(&[inner_size as u64]);

        assert_eq!(inner_generator.pow(&[inner_size as u64]), Fr::one());

        // let precomputed_twiddle_factors = vec![];
        let precomputed_twiddle_factors =  crate::fft::strided_fft::utils::precompute_twiddle_factors(&inner_generator, inner_size as usize);
        non_generic_radix_sqrt::<_, 128>(&mut values_for_sqrt, &omega, &precomputed_twiddle_factors, &worker);

        assert_eq!(values, values_for_sqrt);
    }

    #[test]
    fn test_roots_of_unity(){
        let generator = Fr::root_of_unity();

        let max_power = Fr::S;
        // it means x^(2^max_power) = 1
        println!("max power {}", max_power);
        // let say we need x^4 = x^(2^2) = 1
        // so we need to compute x^(max_power-2)
        let power: u64 = 1<<(max_power - 2);

        let powered_generator = generator.pow(&[power as u64]);

        let actual = powered_generator.pow(&[(1<<2) as u64]);

        assert_eq!(actual, Fr::one());
    }

    #[test]
    fn test_inner_outer_size(){

        let x: usize = 1024;

        // we have 2^10 = 1024 values
        // we want to compute square root size
        // 1. compute log_n by using trailing zeroes
        let log_n = x.trailing_zeros();
        assert_eq!(1<<log_n, x);
        // 2. compute half of log_n
        let half_of_log_n = log_n/2;

        // 3. compupte 2^(logn/2)
        let result = 1<<half_of_log_n;

        assert_eq!(result*result, x);
    }
}
