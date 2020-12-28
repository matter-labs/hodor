use criterion::{Criterion, Bencher};
use rand::{Rng, XorShiftRng, SeedableRng};
use hodor::ff::*;
use criterion_convenience::*;

fn compare_fft_unrolling_16(c: &mut Criterion) {
    use hodor::optimized_fields::naive_f125::Fr;
    let log_2_sizes: Vec<usize> = vec![4, 6, 8, 10, 12, 14, 16, 18, 20, 22];

    type U = (hodor::fft::multicore::Worker, Vec<Fr>, Vec<Fr>, Fr);

    // trivial one, may be change later
    let generator = |log_2_size: &usize| {
        let log_2_size = *log_2_size;
        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let worker = hodor::fft::multicore::Worker::new();
        let size = 1 << log_2_size;
        let values: Vec<Fr> = (0..size).map(|_| rng.gen::<Fr>()).collect();
        let main_domain = hodor::domains::Domain::<Fr>::new_for_size(size as u64).unwrap();
        let omega = main_domain.generator;
        let (inner, outer) = hodor::fft::strided_fft::non_generic::calcualate_inner_and_outer_sizes(size);
        let outer_omega = omega.pow(&[inner as u64]);
        let twiddles = hodor::fft::strided_fft::utils::precompute_twiddle_factors(&outer_omega, outer);

        (worker, values, twiddles, omega)
    };
    let executors = vec![
        ("Unroll 128".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (worker, values, twiddles, omega) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt::<_, 128>(&mut values, &omega, &twiddles, &worker))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        ),
        ("Unroll 256".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (worker, values, twiddles, omega) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt::<_, 256>(&mut values, &omega, &twiddles, &worker))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        ),
        ("Unroll 1024".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (worker, values, twiddles, omega) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt::<_, 1024>(&mut values, &omega, &twiddles, &worker))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        )
    ];

    log_parametrized_comparison_benchmark(c, 
        "FFT for 125 bit field and various unrolls", 
        true,
        &log_2_sizes, 
        generator,
        executors
    );
}

fn compare_fft_unrolling_32(c: &mut Criterion) {
    use hodor::optimized_fields::f252::Fr;
    let log_2_sizes: Vec<usize> = vec![4, 6, 8, 10, 12, 14, 16, 18, 20, 22];

    type U = (hodor::fft::multicore::Worker, Vec<Fr>, Vec<Fr>, Fr);

    // trivial one, may be change later
    let generator = |log_2_size: &usize| {
        let log_2_size = *log_2_size;
        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let worker = hodor::fft::multicore::Worker::new();
        let size = 1 << log_2_size;
        let values: Vec<Fr> = (0..size).map(|_| rng.gen::<Fr>()).collect();
        let main_domain = hodor::domains::Domain::<Fr>::new_for_size(size as u64).unwrap();
        let omega = main_domain.generator;
        let (inner, outer) = hodor::fft::strided_fft::non_generic::calcualate_inner_and_outer_sizes(size);
        let outer_omega = omega.pow(&[inner as u64]);
        let twiddles = hodor::fft::strided_fft::utils::precompute_twiddle_factors(&outer_omega, outer);

        (worker, values, twiddles, omega)
    };
    let executors = vec![
        ("Unroll 128".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (worker, values, twiddles, omega) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt::<_, 128>(&mut values, &omega, &twiddles, &worker))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        ),
        ("Unroll 256".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (worker, values, twiddles, omega) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt::<_, 256>(&mut values, &omega, &twiddles, &worker))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        ),
        ("Unroll 1024".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (worker, values, twiddles, omega) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt::<_, 1024>(&mut values, &omega, &twiddles, &worker))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        )
    ];

    log_parametrized_comparison_benchmark(c, 
        "FFT for 252 bit field and various unrolls", 
        true,
        &log_2_sizes, 
        generator,
        executors
    );
}

pub fn group(crit: &mut Criterion) {
    compare_fft_unrolling_16(crit);
    compare_fft_unrolling_32(crit);
}