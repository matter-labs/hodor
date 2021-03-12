use core::num;
use std::convert::TryInto;

use criterion::{
    measurement::{Measurement, WallTime},
    AxisScale, Bencher, BenchmarkGroup, BenchmarkId, Criterion, PlotConfiguration, Throughput,
};
use hodor::optimized_fields::naive_f125 as Fr125Naive;
use hodor::{fft::fft::best_fft, optimized_fields::f125::Fr as Fr125Asm};
use hodor::{fft::strided_fft::non_generic, optimized_fields::f252::Fr as FrAsm};
use hodor::{
    fft::strided_fft::non_generic::non_generic_radix_sqrt_with_rayon,
    optimized_fields::f252_asm_generated::Fr as FrAsmG,
};

use hodor::{
    domains::Domain,
    ff::{Field, PrimeField},
    fft::{
        multicore::Worker,
        strided_fft::{
            non_generic::calcualate_inner_and_outer_sizes, utils::precompute_twiddle_factors,
        },
    },
};
use hodor::{
    fft::strided_fft::shuffle::{
        transpose_square, transpose_square_with_chunks, transpose_square_with_square_tiles,
    },
    optimized_fields::naive_f252::Fr as FrNaive,
};
use openzkp_primefield::{
    fft::fft_vec_recursive_with_unroll, fft::get_twiddles, fft::radix_sqrt_with_unroll,
    fft::transpose_square_1 as transpose_square_openzkp, u256::U256, Fft, FieldElement, FieldLike,
    Root,
};

use hodor::fft::strided_fft::non_generic::{radix_2_butterfly, radix_2_butterfly_with_twiddle};
use openzkp_primefield::fft::small::{radix_2, radix_2_twiddle};

use hodor::fft::strided_fft::non_generic::non_generic_radix_sqrt;
use rand::{Rand, Rng, SeedableRng, XorShiftRng};

use rayon::{prelude::*, ThreadPoolBuilder};

struct OpenZkpBencher;

impl OpenZkpBencher {
    fn compute_omega() -> FieldElement {
        FieldElement::from(1 as u64)
    }

    fn generate_fft_values(size: usize) -> Vec<FieldElement> {
        let mut values = vec![];

        for i in 0..size {
            let el = FieldElement::from(i);
            values.push(el);
        }

        values
    }

    fn generate_random_values<R: Rng>(size: usize, rng: &mut R) -> Vec<FieldElement> {
        let mut values = vec![];

        for _ in 0..size {
            let el = FieldElement::from(U256::from_limbs(rng.gen()));
            values.push(el);
        }

        values
    }

    fn random_fr<R: Rng>(rng: &mut R) -> FieldElement {
        FieldElement::from(U256::from_limbs(rng.gen()))
    }

    fn bench_fr_mul<R: Rng>(group: &mut BenchmarkGroup<WallTime>, rng: &mut R) {
        let el = Self::random_fr(rng);
        group.bench_with_input("openzkp", &el, |b, elt| b.iter(|| elt * elt));
    }

    // TODO: do we need size?
    fn bench_fft_butterfly<R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        rng: &mut R,
        size: usize,
        offset: usize,
    ) {
        let mut values = Self::generate_random_values(size, rng);

        group.bench_function("openzkp-radix-2", |b| {
            b.iter(|| {
                radix_2(&mut values, offset, 1);
            })
        });
    }

    fn bench_fft_butterfly_with_twiddles<R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        rng: &mut R,
        size: usize,
        offset: usize,
    ) {
        let mut values = Self::generate_random_values(size, rng);
        let twiddle = Self::random_fr(rng);

        group.bench_function("openzkp-radix-2-with-twiddles", |b| {
            b.iter(|| {
                radix_2_twiddle(&mut values, &twiddle, offset, 1);
            })
        });
    }

    fn bench_matrix_transposition<R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        dim: usize,
        rng: &mut R,
    ) {
        let mut matrix = Self::generate_random_values(dim * dim, rng);

        group.bench_function(format!("openzkp-with-dim-{}", dim), |b| {
            b.iter(|| {
                transpose_square_openzkp(&mut matrix, dim);
            })
        });
    }

    // Since we only interested with single type, use concrete Walltime here.
    // fn benchmark_radix_sqrt<M: Measurement, const LOOP_UNROLL: usize>(group: &mut BenchmarkGroup<WallTime>, size: usize){
    fn bench_fft_radix_sqrt<const LOOP_UNROLL: usize>(
        group: &mut BenchmarkGroup<WallTime>,
        size: usize,
    ) {
        let mut values = Self::generate_fft_values(size);
        let omega = Self::compute_omega();
        let id_params = format!("size-{}", size);
        let id_str = format!("openzkp-radix-sqrt-loop-unroll-{}", LOOP_UNROLL);
        let id = BenchmarkId::new(&id_str, &id_params);
        group.bench_function(id, |b| {
            b.iter(|| radix_sqrt_with_unroll::<_, LOOP_UNROLL>(&mut values, &omega))
        });
    }
}

struct MatterBencher;

impl MatterBencher {
    fn generate_values<F: PrimeField>(size: usize) -> (F, Vec<F>, Vec<F>) {
        let domain: Domain<F> = Domain::new_for_size(size as u64).unwrap();
        let (inner_size, outer_size) = calcualate_inner_and_outer_sizes(size);
        let omega = domain.generator;
        let omega_outer = omega.pow(&[inner_size as u64]);
        let twiddles = precompute_twiddle_factors(&omega_outer, outer_size);

        let mut values = vec![];
        for i in 0..size {
            let el = F::from_str(&i.to_string()).unwrap();
            values.push(el);
        }

        (omega, twiddles, values)
    }

    fn generate_random_fr<F: PrimeField, R: Rng>(rng: &mut R) -> F {
        F::rand(rng)
    }

    fn generate_random_values<F: PrimeField, R: Rng>(size: usize, rng: &mut R) -> Vec<F> {
        let mut values = vec![];
        for _ in 0..size {
            let el = F::rand(rng);
            values.push(el);
        }

        values
    }

    fn bench_fr_mul<F: PrimeField, R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        rng: &mut R,
        fr_name: &str,
    ) {
        let mut tmp: F = Self::generate_random_fr(rng);
        let el = Self::generate_random_fr(rng);
        let bench_id = format!("matter-{}", fr_name);
        group.bench_with_input(bench_id, &el, |b, elt| {
            b.iter(|| {
                tmp.mul_assign(elt);
            })
        });
    }

    fn bench_fft_butterfly<F: PrimeField, R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        rng: &mut R,
        size: usize,
        offset: usize,
    ) {
        let mut values: Vec<F> = Self::generate_random_values(size, rng);
        let twiddle: F = Self::generate_random_fr(rng);

        group.bench_function("matter-radix-2", |b| {
            b.iter(|| {
                radix_2_butterfly(&mut values, offset, 1);
            })
        });
    }

    fn bench_fft_butterfly_with_twiddles<F: PrimeField, R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        rng: &mut R,
        size: usize,
        offset: usize,
    ) {
        let mut values = Self::generate_random_values(size, rng);
        let twiddle: F = Self::generate_random_fr(rng);
        let stride = 1; // simply assume 1;

        group.bench_function("matter-radix-2-with-twiddles", |b| {
            b.iter(|| {
                radix_2_butterfly_with_twiddle(&mut values, &twiddle, offset, stride);
            })
        });
    }

    fn bench_matrix_transposition<F: PrimeField, R: Rng>(
        group: &mut BenchmarkGroup<WallTime>,
        dim: usize,
        rng: &mut R,
    ) {
        let size = dim * dim;
        let mut matrix = Self::generate_random_values::<F, _>(size, rng);

        group.bench_function(format!("matter-dim-{}", dim), |b| {
            b.iter(|| {
                transpose_square(&mut matrix, dim);
            })
        });

        group.bench_function(format!("matter-with-chunks-2x2-dim-{}", dim), |b| {
            b.iter(|| transpose_square_with_chunks::<_, 2>(&mut matrix, dim))
        });

        // group.bench_function(format!("matter-with-chunks-4x4-dim-{}", dim), |b| {
        //     b.iter(|| transpose_square_with_chunks::<_, 4>(&mut matrix, dim))
        // });

        group.bench_function(format!("matter-with-tiles-2x2-dim-{}", dim), |b| {
            b.iter(|| transpose_square_with_square_tiles::<_, 2>(&mut matrix, dim))
        });

        // group.bench_function(format!("matter-with-tiles-4x4-dim-{}", dim), |b| {
        //     b.iter(|| transpose_square_with_square_tiles::<_, 4>(&mut matrix, dim))
        // });
    }

    fn bench_non_generic_radix_sqrt_fft<F: PrimeField, const LOOP_UNROLL: usize>(
        group: &mut BenchmarkGroup<WallTime>,
        size: usize,
        worker: &Worker,
    ) {
        let bit_len = F::CAPACITY;

        let id_params = format!("size-{}", size);
        let id_str = format!(
            "matter-radix-sqrt-Fr-{}-loop-unroll-{}",
            bit_len, LOOP_UNROLL
        );
        let id = BenchmarkId::new(&id_str, &id_params);
        let (omega, twiddles, mut values_m): (F, Vec<F>, Vec<F>) = Self::generate_values(size);
        group.bench_function(id, |b| {
            b.iter(|| {
                non_generic_radix_sqrt::<_, LOOP_UNROLL>(&mut values_m, &omega, &twiddles, &worker)
            });
        });
    }

    fn bench_non_generic_radix_sqrt_fft_with_scalar_fields<
        F: PrimeField,
        const LOOP_UNROLL: usize,
    >(
        group: &mut BenchmarkGroup<WallTime>,
        size: usize,
        worker: &Worker,
        scalar_field: &str,
    ) {
        let id_params = format!("size-{}-loop-unroll-{}", size, LOOP_UNROLL);
        let id_str = format!("matter-radix-sqrt-field-{}", scalar_field);
        let id = BenchmarkId::new(&id_str, &id_params);
        let (omega, twiddles, mut values_m): (F, Vec<F>, Vec<F>) = Self::generate_values(size);
        group.bench_function(id, |b| {
            b.iter(|| {
                non_generic_radix_sqrt::<_, LOOP_UNROLL>(&mut values_m, &omega, &twiddles, &worker)
            });
        });
    }
}

fn bench_init(crit: &mut Criterion, group_id: String) -> (BenchmarkGroup<WallTime>, XorShiftRng) {
    let mut group = crit.benchmark_group(group_id);
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    let rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    (group, rng)
}

fn compare_mac(crit: &mut Criterion) {
    use ff::mac_with_carry;
    use openzkp_primefield::u256::mac;

    let (mut group, mut rng) = bench_init(crit, "repr mac".to_string());

    let x: u64 = rng.gen();
    let y: u64 = rng.gen();
    let mut carry = 0u64;

    fn mac_new(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
        let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);

        // #[allow(clippy::cast_possible_truncation)]
        (ret as u64, (ret >> 64) as u64)
    }

    group.bench_function("openzkp mac", |b| {
        b.iter(|| {
            mac(0, x, y, carry);
        })
    });
    group.bench_function("matter mac", |b| {
        b.iter(|| {
            mac_with_carry(0, x, y, &mut carry);
        })
    });
    group.bench_function("matter mac_new", |b| {
        b.iter(|| {
            mac_new(0, x, y, carry);
        })
    });
}

fn compare_ef_fr_multiplication(crit: &mut Criterion) {
    use hodor::optimized_fields::f125::Fr as FrAsm;
    use hodor::optimized_fields::naive_f125::Fr as FrNaive;

    let (mut group, mut rng) = bench_init(crit, "Fr Multiplication(EF)".to_string());

    MatterBencher::bench_fr_mul::<FrAsm, _>(&mut group, &mut rng, "fr-asm");
    MatterBencher::bench_fr_mul::<FrNaive, _>(&mut group, &mut rng, "fr-naive");
    // MatterBencher::bench_fr_mul::<FrAsmG, _>(&mut group, &mut rng, "fr-asm-generated");
    OpenZkpBencher::bench_fr_mul(&mut group, &mut rng);

    group.finish();
}

fn compare_proth_fr_multiplication(crit: &mut Criterion) {
    let (mut group, mut rng) = bench_init(crit, "Fr Multiplication(Proth)".to_string());

    MatterBencher::bench_fr_mul::<FrAsm, _>(&mut group, &mut rng, "fr-asm");
    MatterBencher::bench_fr_mul::<FrNaive, _>(&mut group, &mut rng, "fr-naive");
    // MatterBencher::bench_fr_mul::<FrAsmG, _>(&mut group, &mut rng, "fr-asm-generated");
    OpenZkpBencher::bench_fr_mul(&mut group, &mut rng);

    group.finish();
}

fn compare_fft_butterfly(crit: &mut Criterion) {
    let (mut group, mut rng) = bench_init(crit, "FFT Butterfly".to_string());

    let size = 16;
    let offset: usize = rng.gen::<usize>() % size;

    MatterBencher::bench_fft_butterfly::<FrNaive, _>(
        &mut group,
        &mut rng,
        size.clone(),
        offset.clone(),
    );
    OpenZkpBencher::bench_fft_butterfly(&mut group, &mut rng, size, offset);

    group.finish();
}

fn compare_fft_butterfly_with_twiddles(crit: &mut Criterion) {
    let (mut group, mut rng) = bench_init(crit, "FFT Butterfly with Twiddles".to_string());

    let size = 16;
    let offset: usize = rng.gen::<usize>() % size;

    MatterBencher::bench_fft_butterfly_with_twiddles::<FrNaive, _>(
        &mut group,
        &mut rng,
        size.clone(),
        offset.clone(),
    );
    OpenZkpBencher::bench_fft_butterfly_with_twiddles(&mut group, &mut rng, size, offset);

    group.finish();
}

fn compare_matrix_transposition(crit: &mut Criterion, range: &std::ops::Range<usize>) {
    let (mut group, mut rng) = bench_init(crit, "nxn Matrix Transposition".to_string());

    for size in range.clone().step_by(2) {
        let dim = 1 << (size / 2);
        group.throughput(criterion::Throughput::Elements(size as u64));

        MatterBencher::bench_matrix_transposition::<FrNaive, _>(&mut group, dim, &mut rng);
        OpenZkpBencher::bench_matrix_transposition(&mut group, dim, &mut rng);
    }
}

fn compare_fft_by_unroll_params(crit: &mut Criterion, range: &std::ops::Range<usize>) {
    let worker = Worker::new();
    let (mut group, _) = bench_init(crit, "FFT".to_string());

    let log_sizes: Vec<usize> = range.to_owned().step_by(2).collect();
    let sizes: Vec<usize> = log_sizes.iter().map(|el| 1 << el).collect();

    for (size, _log_size) in sizes.iter().cloned().zip(log_sizes) {
        group.throughput(criterion::Throughput::Elements(size as u64));

        MatterBencher::bench_non_generic_radix_sqrt_fft::<FrAsm, 128>(&mut group, size, &worker);

        MatterBencher::bench_non_generic_radix_sqrt_fft::<FrAsm, 1024>(&mut group, size, &worker);

        OpenZkpBencher::bench_fft_radix_sqrt::<128>(&mut group, size);

        OpenZkpBencher::bench_fft_radix_sqrt::<1024>(&mut group, size);
    }
}

fn compare_fft_by_fields_different_bits(crit: &mut Criterion, range: &std::ops::Range<usize>) {
    let worker = Worker::new();
    let (mut group, _) = bench_init(crit, "FFT".to_string());

    let log_sizes: Vec<usize> = range.to_owned().step_by(2).collect();
    let sizes: Vec<usize> = log_sizes.iter().map(|el| 1 << el).collect();

    for (size, _log_size) in sizes.iter().cloned().zip(log_sizes) {
        group.throughput(criterion::Throughput::Elements(size as u64));

        MatterBencher::bench_non_generic_radix_sqrt_fft::<FrAsm, 128>(&mut group, size, &worker);

        MatterBencher::bench_non_generic_radix_sqrt_fft::<Fr125Asm, 128>(&mut group, size, &worker);
    }
}

fn compare_fft_by_scalar_fields(crit: &mut Criterion, range: &std::ops::Range<usize>) {
    let worker = Worker::new();
    let (mut group, _) = bench_init(crit, "FFT by Different Fr".to_string());

    let log_sizes: Vec<usize> = range.to_owned().step_by(2).collect();
    let sizes: Vec<usize> = log_sizes.iter().map(|el| 1 << el).collect();

    for (size, _log_size) in sizes.iter().cloned().zip(log_sizes) {
        group.throughput(criterion::Throughput::Elements(size as u64));

        MatterBencher::bench_non_generic_radix_sqrt_fft_with_scalar_fields::<FrNaive, 128>(
            &mut group, size, &worker, "fr-naive",
        );

        MatterBencher::bench_non_generic_radix_sqrt_fft_with_scalar_fields::<FrAsm, 128>(
            &mut group, size, &worker, "fr-asm",
        );

        MatterBencher::bench_non_generic_radix_sqrt_fft_with_scalar_fields::<FrAsm, 128>(
            &mut group,
            size,
            &worker,
            "fr-default",
        );

        OpenZkpBencher::bench_fft_radix_sqrt::<128>(&mut group, size);
    }
}

fn compare_square_root_fft_with_best_fft(crit: &mut Criterion, range: std::ops::Range<usize>) {
    let worker = Worker::new();
    let mut group = crit.benchmark_group("FFT comparison");
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    let field_bit_size = FrAsm::CAPACITY;
    println!("running for field {}", field_bit_size);

    for log_size in range {
        let size = 1 << log_size;

        const LOOP_UNROLL: usize = 128;

        let id_params = format!("size-{}", size);
        let id = BenchmarkId::new("radix-sqrt-field-unroll-128", &id_params);
        let (omega, twiddles, mut values): (FrAsm, Vec<FrAsm>, Vec<FrAsm>) =
            MatterBencher::generate_values(size);
        println!("size of elements {}", size);

        let mut values_for_square_root = values.clone();
        group.bench_function(id, |b| {
            b.iter(|| {
                non_generic_radix_sqrt::<_, LOOP_UNROLL>(
                    &mut values_for_square_root,
                    &omega,
                    &twiddles,
                    &worker,
                )
            });
        });

        let id_params = format!("size-{}", size);
        let id = BenchmarkId::new("best-fft", &id_params);
        group.bench_function(id, |b| {
            b.iter(|| {
                best_fft(&mut values, &worker, &omega, log_size as u32, None);
            });
        });
    }
}

fn bench_non_generic_radix_sqrt_fft_with_different_cpus(crit: &mut Criterion) {
    let mut group = crit.benchmark_group("Square root fft with CPUs");
    for num_cpu in (8..50).step_by(2) {
        MatterBencher::bench_non_generic_radix_sqrt_fft::<FrAsm, 128>(
            &mut group,
            1 << 22,
            &Worker::new_with_cpus(num_cpu),
        );
    }
}

fn bench_fft_with_custom_num_threads(crit: &mut Criterion, log_size: usize) {
    let mut group = crit.benchmark_group(format!("FFT(2^{}) Comparison by Num CPUs", log_size));
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Linear);
    // let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    let size = 1 << log_size;

    let mut values_o = OpenZkpBencher::generate_fft_values(size);
    let omega_o = OpenZkpBencher::compute_omega();

    type F = FrAsm;
    let (omega_m, twiddles, mut values_m): (F, Vec<F>, Vec<F>) =
        MatterBencher::generate_values(size);
    for num_threads in vec![2, 4, 8, 16, 24, 32, 48] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("some pool");

        group.bench_function(
            BenchmarkId::new(format!("openzkp-size-{}-unroll-128", size), num_threads),
            |b| {
                b.iter(|| {
                    pool.install(|| radix_sqrt_with_unroll::<_, 128>(&mut values_o, &omega_o));
                })
            },
        );

        let worker = Worker::new_with_cpus(num_threads);
        group.bench_function(
            BenchmarkId::new(format!("matter-size-{}-unroll-128", size), num_threads),
            |b| {
                b.iter(|| {
                    non_generic_radix_sqrt::<_, 128>(&mut values_m, &omega_m, &twiddles, &worker)
                });
            },
        );

        group.bench_function(
            BenchmarkId::new(format!("best-fft-{}", size), num_threads),
            |b| {
                b.iter(|| {
                    best_fft(&mut values_m, &worker, &omega_m, log_size as u32, None);
                });
            },
        );
    }
}

fn compare_matrix_transposition_with_custom_num_threads(crit: &mut Criterion) {
    use rayon::ThreadPoolBuilder;

    let mut group = crit.benchmark_group("Matrix Transposition Comparison by Num CPUs");
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Linear);
    group.plot_config(new_config);

    let log_sizes = vec![4, 8];
    // let log_sizes = vec![22, 24, 26];
    for log_size in log_sizes {
        let dim = 1 << (log_size / 2);

        let size = 1 << log_size;

        let mut values_o = OpenZkpBencher::generate_fft_values(size);
        let omega_o = OpenZkpBencher::compute_omega();

        type F = FrAsm;
        let (omega_m, twiddles, mut values_m): (F, Vec<F>, Vec<F>) =
            MatterBencher::generate_values(size);
        // for num_threads in vec![2, 4, 8, 16, 24, 32, 48] {
        for num_threads in vec![2, 4] {
            let pool = ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .expect("some pool");

            group.bench_function(
                BenchmarkId::new(format!("openzkp-size-{}", size), num_threads),
                |b| {
                    b.iter(|| {
                        pool.install(|| transpose_square_openzkp(&mut values_o, dim));
                    })
                },
            );

            // let worker = Worker::new_with_cpus(num_threads);
            group.bench_function(
                BenchmarkId::new(format!("matter-size-{}", size), num_threads),
                |b| {
                    b.iter(|| {
                        transpose_square(&mut values_m, dim);
                    });
                },
            );

            group.bench_function(
                BenchmarkId::new(format!("matter-with-tiles-size-{}", size), num_threads),
                |b| {
                    b.iter(|| {
                        transpose_square_with_square_tiles::<_, 2>(&mut values_m, dim);
                    });
                },
            );

            group.bench_function(
                BenchmarkId::new(format!("matter-with-chunks-fft-{}", size), num_threads),
                |b| {
                    b.iter(|| {
                        transpose_square_with_chunks::<_, 2>(&mut values_m, dim);
                    });
                },
            );
        }
    }
}

fn bench_simultaneous_fft_with_num_threads(
    crit: &mut Criterion,
    fft_log_size: usize,
    number_of_simulataneous_fft: usize,
) {
    let size = 1 << fft_log_size;

    let mut group = crit.benchmark_group(format!(
        "{} Simultaneous FFT(1<<{}) by CPUs",
        number_of_simulataneous_fft, fft_log_size
    ));

    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    type F = FrAsm;
    let (omega_m, twiddles, values_m): (F, Vec<F>, Vec<F>) = MatterBencher::generate_values(size);
    let worker = Worker::new();

    // we need to run for 1,2,3,4,5,6 parralel jobs fft

    let mut values = values_m.clone();
    for _ in 1..number_of_simulataneous_fft {
        values.extend_from_slice(&values_m);
    }

    // for num_cpus in vec![6,24,48] {
    for num_cpus in vec![4, 8, 16, 24, 32, 48] {
        if num_cpus < number_of_simulataneous_fft || num_cpus % number_of_simulataneous_fft != 0 {
            continue;
        }
        // group.bench_function(
        //     BenchmarkId::new("sqrt-simultaneous-fft", format!("{}-cpu", num_cpus)),
        //     |b| {
        //         b.iter(|| {
        //             run_simultaneous_fft(
        //                 &mut values,
        //                 &omega_m,
        //                 &twiddles,
        //                 &worker,
        //                 num_cpus,
        //                 number_of_simulataneous_fft,
        //             );
        //         });
        //     },
        // );
    }
}

fn bench_simultaneous_fft_with_worker(crit: &mut Criterion, log2_size: usize) {
    let size = 1 << log2_size;

    let mut group = crit.benchmark_group(format!("[Matter-with-Worker]Simultaneous FFTs(1<<{}) by CPUs", log2_size));

    type F = FrAsm;
    let (omega, twiddles, values_m): (F, Vec<F>, Vec<F>) = MatterBencher::generate_values(size);

    for num_cpus in vec![4, 8, 12, 16, 20, 24, 32, 36, 40, 48] {
        
        for num_fft in vec![1, 2, 4, 8, 12, 16, 24, 32, 48] {
            if num_fft > num_cpus {
                continue;
            }
            let mut values = values_m.clone();
            for _ in 1..num_fft {
                values.extend_from_slice(&values_m);
            }
            let bench_id = BenchmarkId::new(format!("{}-simultaneous-fft", num_fft), num_cpus);
            group.bench_function(bench_id, |b| {
                b.iter(|| run_simultaneous_fft(&mut values, &omega, &twiddles, num_cpus, num_fft));
            });
        }
    }
}

fn run_simultaneous_fft<F: PrimeField>(
    values: &mut [F],
    omega: &F,
    twiddles: &[F],
    total_num_cpus: usize,
    number_of_simulataneous_fft: usize,
) {
    let worker = Worker::new_with_cpus(total_num_cpus);

    let cpu_per_fft = total_num_cpus / number_of_simulataneous_fft;

    let number_of_values_per_fft = values.len() / number_of_simulataneous_fft;

    let mut values_ref = values.as_mut();

    worker.scope(0, |scope, _| {
        for _ in 0..number_of_simulataneous_fft {
            let inner_worker = Worker::new_with_cpus(cpu_per_fft);
            let (values_for_thread, rest) = values_ref.split_at_mut(number_of_values_per_fft);
            assert!(values_for_thread.len().is_power_of_two());
            values_ref = rest;
            scope.spawn(move |_| {
                assert_eq!(values_for_thread.len(), number_of_values_per_fft);
                non_generic_radix_sqrt::<_, 128>(
                    values_for_thread,
                    &omega,
                    &twiddles,
                    &inner_worker,
                );
            });
        }
    });
}

fn run_simultaneous_fft_for_matter_with_rayon<F: PrimeField>(
    values: &mut [F],
    omega: &F,
    twiddles: &[F],
    total_num_cpus: usize,
    number_of_simulataneous_fft: usize,
) {
    let number_of_values_per_fft = values.len() / number_of_simulataneous_fft;

    values
        .par_chunks_mut(number_of_values_per_fft)
        .for_each(|child_values| {
            non_generic_radix_sqrt_with_rayon::<_, 128>(child_values, omega, twiddles);
        });
}

use num_traits::{Inv, Pow};
use openzkp_primefield::FieldOps;

fn run_simultaneous_fft_for_openzkp<F>(
    values: &mut [F],
    omega: &F,
    total_num_cpus: usize,
    number_of_simulataneous_fft: usize,
) where
    F: FieldLike + Send + Sync,
    for<'a, 'r> &'a F:
        Pow<usize, Output = F> + Inv<Output = Option<F>> + FieldOps<F, F> + FieldOps<&'r F, F>,
{
    let number_of_values_per_fft = values.len() / number_of_simulataneous_fft;

    values
        .par_chunks_mut(number_of_values_per_fft)
        .for_each(|child_values| {
            radix_sqrt_with_unroll::<_, 128>(child_values, omega);
        });
}

fn bench_simultaneous_fft_for_openzkp(crit: &mut Criterion, log2_size: usize) {
    let size = 1 << log2_size;

    let mut group = crit.benchmark_group(format!(
        "[OpenZKP] Simultaneous FFT(1<<{}) by CPUS ",
        log2_size
    ));

    let values_o = OpenZkpBencher::generate_fft_values(size);
    let omega = OpenZkpBencher::compute_omega();

    for num_cpus in vec![4, 8, 12, 16, 20, 24, 32, 36, 40, 48] {
        // for num_cpus in vec![4, 8] {
        // for num_fft in vec![2,4] {
        for num_fft in vec![1, 2, 4, 8, 12, 16, 24, 32, 48] {
            if num_fft > num_cpus {
                continue;
            }
            let mut values = values_o.clone();
            for _ in 1..num_fft {
                values.extend_from_slice(&values_o);
            }

            let bench_id = BenchmarkId::new(format!("{}-simultaneous-fft", num_fft), num_cpus);
            group.bench_function(bench_id, |b| {
                let pool = ThreadPoolBuilder::new()
                    .num_threads(num_cpus)
                    .build()
                    .expect("some pool");

                b.iter(|| {
                    pool.install(|| {
                        run_simultaneous_fft_for_openzkp(&mut values, &omega, num_cpus, num_fft)
                    })
                });
            });
        }
    }
}

fn bench_simultaneous_fft_for_matter_with_rayon(crit: &mut Criterion, log2_size: usize) {
    let size = 1 << log2_size;

    let mut group = crit.benchmark_group(format!(
        "[Matter] Simultaneous FFT(1<<{}) by CPUS ",
        log2_size
    ));

    type F = FrAsm;
    let (omega, twiddles, values_m): (F, Vec<F>, Vec<F>) = MatterBencher::generate_values(size);

    for num_cpus in vec![4, 8, 12, 16, 20, 24, 32, 36, 40, 48] {

        for num_fft in vec![1, 2, 4, 8, 12, 16, 24, 32, 48] {
            if num_fft > num_cpus {
                continue;
            }
            let mut values = values_m.clone();
            for _ in 1..num_fft {
                values.extend_from_slice(&values_m);
            }

            let bench_id = BenchmarkId::new(format!("{}-simultaneous-fft", num_fft), num_cpus);
            group.bench_function(bench_id, |b| {
                let pool = ThreadPoolBuilder::new()
                    .num_threads(num_cpus)
                    .build()
                    .expect("some pool");

                b.iter(|| {
                    pool.install(|| {
                        run_simultaneous_fft_for_matter_with_rayon(
                            &mut values,
                            &omega,
                            &twiddles,
                            num_cpus,
                            num_fft,
                        )
                    })
                });
            });
        }
    }
}

pub fn group(crit: &mut Criterion) {
    // compare_mac(crit);
    // compare_proth_fr_multiplication(crit);
    // compare_ef_fr_multiplication(crit);
    //     compare_fft_butterfly(crit);
    //     compare_fft_butterfly_with_twiddles(crit);
    // let matrix_range = 16..24;
    // compare_matrix_transposition(crit, &matrix_range);
    // for log_size in vec![22, 24, 26]{
    //     bench_fft_with_custom_num_threads(crit, log_size);
    // }
    bench_simultaneous_fft_with_worker(crit, 22);
    for log2_size in vec![24, 26] {
        bench_simultaneous_fft_with_worker(crit, log2_size);
        bench_simultaneous_fft_for_openzkp(crit, log2_size);
        bench_simultaneous_fft_for_matter_with_rayon(crit, log2_size);
    }
    // bench_simultaneous_fft_with_worker(crit, 22);
    // bench_simultaneous_fft_for_openzkp(crit, 22);
    // bench_simultaneous_fft_for_matter_with_rayon(crit, 22);
    // compare_matrix_transposition_with_custom_num_threads(crit);
    // bench_two_simultaneous_fft_with_custom_num_threads(crit);
    // let fft_range = 8..12;
    // compare_fft_by_unroll_params(crit, &fft_range);
    // compare_fft_by_fields_different_bits(crit, &fft_range);
    //     compare_fft_by_scalar_fields(crit, &fft_range);
    // compare_square_root_fft_with_best_fft(crit, fft_range);
}
