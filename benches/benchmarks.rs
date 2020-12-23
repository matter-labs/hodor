extern crate ff;
extern crate rand;
extern crate hodor;

use self::ff::*;

mod fr {
    use crate::ff::*;

    #[derive(PrimeField)]
    #[PrimeFieldModulus = "63802944035360449460622495747797942273"]
    #[PrimeFieldGenerator = "3"]
    pub struct Fr(FrRepr);
}

mod fr256 {
    use crate::ff::*;
    #[derive(PrimeField)]
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    pub struct Fr(FrRepr);
}

use criterion::{black_box, criterion_group, criterion_main, Criterion, Bencher, PlotConfiguration, AxisScale, Throughput, BenchmarkId};
use std::convert::*;

fn bench_for_parameter_set<F, T: Copy + Clone + TryInto<u64> + 'static>(c: &mut Criterion, id: &str, with_throughput: bool, parameters: &[T], mut f: F)
where
    F: FnMut(&mut Bencher<'_>, &T),
    <T as TryInto<u64>>::Error: std::fmt::Debug
{
    let mut group = c.benchmark_group(id);
    for param in parameters.iter() {
        let param_as_u64: u64 = param.clone().try_into().unwrap_or_else(|e| panic!(format!("invalid parameter, can not convert to u64: {:?}", e)));
        if with_throughput {
            group.throughput(Throughput::Elements(param_as_u64));
        }
        // let pretty_description = format!("Parameter = {}", param_as_u64);
        group.bench_with_input(BenchmarkId::from_parameter(param_as_u64), param, |b, param| {
            f(b, param)
        });
    }
    group.finish();
}

fn bench_for_log_parameter_set<F, T: Copy + Clone + TryInto<u64> + 'static>(c: &mut Criterion, id: &str, with_throughput: bool, parameter_logs: &[T], mut f: F)
where
    F: FnMut(&mut Bencher<'_>, &T),
    <T as TryInto<u64>>::Error: std::fmt::Debug
{
    let mut group = c.benchmark_group(id);
    let plot_config = PlotConfiguration::default()
        .summary_scale(AxisScale::Logarithmic);
    group.plot_config(plot_config);
    for param in parameter_logs.iter() {
        let log_param_as_u64: u64 = param.clone().try_into().unwrap_or_else(|e| panic!(format!("invalid parameter, can not convert to u64: {:?}", e)));
        let param_as_u64: u64 = 1 << log_param_as_u64;
        if with_throughput {
            group.throughput(Throughput::Elements(param_as_u64));
        }
        // let pretty_description = format!("Parameter = {}, log2 of parametter = {}", param_as_u64, log_param_as_u64);
        let pretty_description = format!("2^{} = {}", log_param_as_u64, param_as_u64);
        group.bench_with_input(BenchmarkId::from_parameter(pretty_description), &param, |b, param| {
            f(b, param)
        });
    }
    group.finish();
}

fn log_parametrized_comparison_benchmark<T: Copy + Clone + TryInto<u64> + 'static, U, P>(
    c: &mut Criterion, 
    id: &str, 
    with_throughput: bool, 
    parameter_logs: &[T], 
    preparation_function: P,
    mut benchmarking_functions: Vec<(String, Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>)>,
) where
    P: Fn(&T) -> U,
    <T as TryInto<u64>>::Error: std::fmt::Debug
{
    let mut group = c.benchmark_group(id);
    let plot_config = PlotConfiguration::default()
        .summary_scale(AxisScale::Logarithmic);
    group.plot_config(plot_config);
    for param in parameter_logs.iter() {
        let log_param_as_u64: u64 = param.clone().try_into().unwrap_or_else(|e| panic!(format!("invalid parameter, can not convert to u64: {:?}", e)));
        let param_as_u64: u64 = 1 << log_param_as_u64;
        if with_throughput {
            group.throughput(Throughput::Elements(param_as_u64));
        }
        let parameters = preparation_function(param);
        for (name, executor) in benchmarking_functions.iter_mut() {
            let pretty_description = format!("2^{} = {}", log_param_as_u64, param_as_u64);
            let id = BenchmarkId::new(&*name, pretty_description);
            group.bench_with_input(id, &parameters, |b, param| {
                executor(b, param)
            });
        }

    }
    group.finish();
}

fn mul_assing_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};
    use fr::Fr;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Mont mul assign 128", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

fn fft_rec_small(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};
    use fr::Fr;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    const SIZE: usize = 1 << 12;
    let domain = hodor::domains::Domain::<Fr>::new_for_size(SIZE as u64).unwrap();
    let twiddles = hodor::fft::strided_fft::utils::precompute_twiddle_factors(&domain.generator, SIZE);
    let values: Vec<Fr> = (0..SIZE).map(|_| rng.gen::<Fr>()).collect();
    let mut values: [Fr; SIZE] = values.try_into().unwrap();
    c.bench_function("Small recursive fft", 
        |bencher| bencher.iter(
            || hodor::fft::strided_fft::fft::small_size_serial_fft::<Fr, SIZE, 128>(&mut values, &twiddles, 0, 1, 1)
        )
    );
}

fn transpose_square_16(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};
    use fr::Fr;

    let rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let log_2_sizes: Vec<usize> = vec![4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26];

    bench_for_log_parameter_set(c, 
        "Transpose square for 16 byte values", 
        true,
        &log_2_sizes, 
        |bencher, &log_2_size| {
            let mut rng = rng.clone();
            let mut values: Vec<Fr> = (0..(1<<log_2_size)).map(|_| rng.gen::<Fr>()).collect();
            let size = 1 << (log_2_size/2);
            bencher.iter(|| hodor::fft::strided_fft::shuffle::transpose_square(&mut values, size))
        }
    );
}

fn transpose_square_32(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};
    use fr256::Fr;

    let rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let log_2_sizes: Vec<usize> = vec![4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26];

    bench_for_log_parameter_set(c, 
        "Transpose square for 32 byte values", 
        true,
        &log_2_sizes, 
        |bencher, &log_2_size| {
        let mut rng = rng.clone();
        let mut values: Vec<Fr> = (0..(1<<log_2_size)).map(|_| rng.gen::<Fr>()).collect();
        let size = 1 << (log_2_size/2);
        bencher.iter(|| hodor::fft::strided_fft::shuffle::transpose_square(&mut values, size))
        }
    );
}

fn compare_fft_unrolling(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};
    use fr::Fr;
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
        "FFT for 128 bit field and various unrolls", 
        true,
        &log_2_sizes, 
        generator,
        executors
    );
}

// criterion_group!(benches, mul_assing_benchmark, fft_rec_small, transpose_square_16, transpose_square_32, fft_sqrt_strategy);
criterion_group!(benches, mul_assing_benchmark, fft_rec_small, compare_fft_unrolling);
criterion_main!(benches);

// criterion_group!(
// 	name = advanced;
//     config = Criterion::default().warm_up_time(std::time::Duration::from_secs(5));
//     targets = mul_assing_benchmark
// );
// criterion_main!(advanced);
