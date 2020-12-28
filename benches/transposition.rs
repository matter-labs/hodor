use rand::{Rng, XorShiftRng, SeedableRng};
use criterion::{Criterion, Bencher};
use criterion_convenience::*;

fn transpose_square_16(c: &mut Criterion) {
    use hodor::optimized_fields::naive_f125::Fr;

    let log_2_sizes: Vec<usize> = vec![4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26];

    type U = (Vec<Fr>, usize);

    // trivial one, may be change later
    let generator = |log_2_size: &usize| {
        let log_2_size = *log_2_size;
        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let size = 1 << log_2_size;
        let values: Vec<Fr> = (0..size).map(|_| rng.gen::<Fr>()).collect();

        let size = 1 << (log_2_size/2);

        (values, size)
    };
    let executors = vec![
        ("Baseline".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (values, size) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::shuffle::transpose_square(&mut values, *size))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        ),
    ];

    log_parametrized_comparison_benchmark(c, 
        "Transpose square for 16 byte values", 
        true,
        &log_2_sizes, 
        generator,
        executors
    );
}

fn transpose_square_32(c: &mut Criterion) {
    use hodor::optimized_fields::f252::Fr;

    let log_2_sizes: Vec<usize> = vec![4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26];

    type U = (Vec<Fr>, usize);

    // trivial one, may be change later
    let generator = |log_2_size: &usize| {
        let log_2_size = *log_2_size;
        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let size = 1 << log_2_size;
        let values: Vec<Fr> = (0..size).map(|_| rng.gen::<Fr>()).collect();

        let size = 1 << (log_2_size/2);

        (values, size)
    };
    let executors = vec![
        ("Baseline".to_string(),
            Box::new(|bencher: &mut Bencher<'_>, params: &U| {
                let (values, size) = &params;
                let mut values = values.to_vec();
                bencher.iter(|| hodor::fft::strided_fft::shuffle::transpose_square(&mut values, *size))
            }) as Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>
        ),
    ];

    log_parametrized_comparison_benchmark(c, 
        "Transpose square for 32 byte values", 
        true,
        &log_2_sizes, 
        generator,
        executors
    );
}

pub fn group(crit: &mut Criterion) {
    transpose_square_16(crit);
    transpose_square_32(crit);
}