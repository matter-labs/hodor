use ff::PrimeField;
use super::utils::*;

pub fn small_size_serial_fft<F: PrimeField, const N: usize, const MAX_LOOP_UNROLL: usize>(
    values: &mut [F; N],
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
            small_size_serial_fft::<F, N, MAX_LOOP_UNROLL>(values, precomputed_twiddle_factors, offset, 2 * count, 2 * stride);
        } else {
            // we may parallelize this too as indexes do not overlap
            small_size_serial_fft::<F, N, MAX_LOOP_UNROLL>(values, precomputed_twiddle_factors, offset, count, 2 * stride);
            small_size_serial_fft::<F, N, MAX_LOOP_UNROLL>(values, precomputed_twiddle_factors, offset + stride, count, 2 * stride);
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