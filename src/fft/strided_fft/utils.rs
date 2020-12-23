use ff::PrimeField;
use super::super::multicore::Worker;

pub const fn bitreverse_index(size: usize, index: usize) -> usize {
    const USIZE_BITS: usize = 0_usize.count_zeros() as usize;
    // debug_assert!(index < size);
    if size == 1 {
        0
    } else {
        // debug_assert!(size.is_power_of_two());
        let bits = size.trailing_zeros() as usize;
        index.reverse_bits() >> (USIZE_BITS - bits)
    }
}

pub fn bitreverse_index_for_constant_size<const N: usize>(index: usize) -> usize {
    const USIZE_BITS: usize = 0_usize.count_zeros() as usize;
    debug_assert!(N.is_power_of_two());
    debug_assert!(index < N);
    if N == 1 {
        0
    } else {
        index.reverse_bits() >> (USIZE_BITS - N)
    }
}

pub fn bitreverse_enumeration<T: Sized>(v: &mut [T]) {
    let n = v.len();
    for i in 0..n {
        let j = bitreverse_index(n, i);
        if j > i {
            v.swap(i, j);
        }
    }
}

#[inline(always)]
pub fn radix_2_butterfly<F: PrimeField, const N: usize>(values: &mut [F; N], offset: usize, stride: usize)
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
pub fn radix_2_butterfly_with_twiddle<F: PrimeField, const N: usize>(values: &mut [F; N], twiddle: &F, offset: usize, stride: usize)
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

pub fn precompute_twiddle_factors<F: PrimeField>(omega: &F, size: usize) -> Vec<F>
{
    debug_assert!(size.is_power_of_two());
    debug_assert_eq!(omega.pow(&[size as u64]), F::one());
    // let mut twiddles = vec![F::zero(); size/2];
    let mut twiddles = (0..size / 2).map(|i| omega.pow(&[i as u64])).collect::<Vec<_>>();
    bitreverse_enumeration(&mut twiddles);
    twiddles
}

// make parallel version
pub fn precompute_twiddle_factors_parallelized<F: PrimeField>(omega: &F, size: usize, worker: &Worker) -> Vec<F>
{
    debug_assert!(size.is_power_of_two());
    debug_assert_eq!(omega.pow(&[size as u64]), F::one());
    let mut twiddles = vec![F::zero(); size/2];
    worker.scope(twiddles.len(), |scope, chunk_size| {
        for (chunk_idx, twiddles) in twiddles.chunks_mut(chunk_size).enumerate() {
            scope.spawn(move |_| {
                let mut power_of_omega = omega.pow(&[(chunk_idx * chunk_size) as u64]);
                for twiddle in twiddles.iter_mut() {
                    *twiddle = power_of_omega;

                    power_of_omega.mul_assign(&omega);
                }
            });
        }
    });

    bitreverse_enumeration(&mut twiddles);
    twiddles
}