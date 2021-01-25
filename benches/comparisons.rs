use criterion::{AxisScale, Bencher, BenchmarkId, Criterion, PlotConfiguration, Throughput};
use criterion_convenience::log_parametrized_comparison_benchmark;
use criterion_convenience::*;
use hodor::{ff, fft::strided_fft::{self, shuffle::{self, transpose_square_with_square_tiles}}, optimized_fields::naive_f125::Fr};
use hodor::{ff::*, fft::multicore, optimized_fields};
use multicore::Worker;
use openzkp_primefield::fft::radix_sqrt;
use openzkp_primefield::fft::transpose_square_stretch;
use openzkp_primefield::{FieldElement, FieldLike};
use rand::{Rng, SeedableRng, XorShiftRng};
use shuffle::{transpose_square_with_chunks, natural_to_morton, recursive_morton, recursive_morton_with_multicore};


fn compare_transpose(crit: &mut Criterion) {
    let sizes = (9..10).map(|e| e * 2).collect::<Vec<usize>>();
    // let sizes = (9..20).collect::<Vec<usize>>();

    let generate_matrix = |log_2_size: &usize| {
        // let log_2_size = *log_2_size;
        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let size = 1 << log_2_size;
        let values: Vec<Fr> = (0..size).map(|_| rng.gen::<Fr>()).collect();
        let dimension_length = 1 << (log_2_size >> 1);
        (values, dimension_length)
    };

    // type BENCH_INPUTS = (Vec<Fr>, usize);

    let mut group = crit.benchmark_group("transpose");
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(plot_config);
    
    let worker = Worker::new();

    for size in sizes.iter() {
        let bench_inputs = generate_matrix(size);


        let title = format!("with-chunks-{}", size);
        let bench_id = BenchmarkId::new("matter", title);
        group.throughput(Throughput::Elements(1 << size));
        group.bench_with_input(bench_id, &bench_inputs, |b, params| {
            let mut matrix = params.0.clone();
            let dim = bench_inputs.1;

            b.iter(|| {
                hodor::fft::strided_fft::shuffle::transpose_square_with_chunks::<_, 4>(
                    &mut matrix,
                    dim,
                )
            })
        });

        let title = format!("with-tiles-4-size-{}", size);
        let bench_id = BenchmarkId::new("matter", title);

        group.bench_with_input(bench_id, &bench_inputs, |b, params| {
            let mut matrix = params.0.clone();
            let dim = bench_inputs.1;
            b.iter(|| transpose_square_with_square_tiles::<_, 4>(&mut matrix, dim))
        });

        let title = format!("with-tiles-8-size-{}", size);
        let bench_id = BenchmarkId::new("matter", title);

        group.bench_with_input(bench_id, &bench_inputs, |b, params| {
            let mut matrix = params.0.clone();
            let dim = bench_inputs.1;
            b.iter(|| transpose_square_with_square_tiles::<_, 8>(&mut matrix, dim))
        });

        let title = format!("with-morton-{}", size);
        let bench_id = BenchmarkId::new("matter", title);

        group.bench_with_input(bench_id, &bench_inputs, |b, params| {
            let mut matrix = params.0.clone();
            let dim = bench_inputs.1;
            b.iter(|| recursive_morton(&mut matrix, dim))
        });
        let title = format!("with-morton-multicore-{}", size);
        let bench_id = BenchmarkId::new("matter", title);

        group.bench_with_input(bench_id, &bench_inputs, |b, params| {
            let mut matrix = params.0.clone();
            let dim = bench_inputs.1;
            b.iter(|| recursive_morton_with_multicore(&mut matrix, dim, &worker))
        });

        let title = format!("open-zkp-{}", size);
        let bench_id = BenchmarkId::new("openzkp", title);

        group.bench_with_input(bench_id, &bench_inputs, |b, params| {
            let mut matrix = params.0.clone();
            let dim = bench_inputs.1;
            // TODO: openzkp fails whe dim = 1<<19, 1 << 17 odd values
            b.iter(|| transpose_square_stretch(&mut matrix, dim, 1))
        });
    }

    group.finish();
}

fn bench_natural_to_morton(crit: &mut Criterion){
    let log_size = 20;
    let dim = 1 << (10);

    let mut matrix = vec![];
    
    for i in 0..(dim * dim) {
        let el = Fr::from_str(&i.to_string()).unwrap();
        matrix.push(el);
    }
    crit.bench_function("natural to morton", |b| b.iter(|| natural_to_morton(&mut matrix, dim)));
}

fn bench_recursive_morton(crit: &mut Criterion){
    let log_size = 20;
    let dim = 1 << (log_size/2);

    let mut matrix = vec![];
    
    for i in 0..(dim * dim) {
        let el = Fr::from_str(&i.to_string()).unwrap();
        matrix.push(el);
    }
    crit.bench_function("natural to morton", |b| b.iter(|| recursive_morton(&mut matrix, dim)));
}

pub fn group(crit: &mut Criterion) {
    compare_transpose(crit);
}
