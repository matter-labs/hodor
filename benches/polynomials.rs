use criterion::{AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use hodor::fft::multicore::Worker;
use hodor::{domains::Domain, optimized_fields::f252_asm_generated::Fr, polynomials::Polynomial};
use hodor::{
    ff::{Field, PrimeField},
    fft::cooley_tukey_ntt::{BitReversedOmegas, CTPrecomputations},
};
use rand::{Rand, Rng, SeedableRng, XorShiftRng};

pub fn group(crit: &mut Criterion) {
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    // let worker = Worker::new();

    let mut group = crit.benchmark_group("LDE");
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    // for num_cpus in vec![1, 2, 4, 8, 16, 24, 32, 48] {
    for num_cpus in vec![2,3,4,5,6,7, 8] {
        let worker = Worker::new_with_cpus(num_cpus);
        let log_size = 10;
        let size: usize = 1 << log_size;
        let lde_factor = 2;
        let coeffs: Vec<_> = (0..size).map(|_| Fr::rand(rng)).collect();
        let poly1 = Polynomial::from_coeffs(coeffs).expect("some poly");
        let poly2 = poly1.clone();

        let omegas_bitreversed =
            BitReversedOmegas::<Fr>::new_for_domain_size(size.next_power_of_two());

        let coset_factor = Fr::one();

        let bench_id = BenchmarkId::new(format!("lde-with-best-fft-{}", size), num_cpus);
        group.bench_function(bench_id, |b| {
            b.iter(|| {
                // poly1.clone().lde(&worker, lde_factor).expect("some lde");
                poly1
                    .clone()
                    .bitreversed_lde_using_bitreversed_ntt(
                        &worker,
                        lde_factor,
                        &omegas_bitreversed,
                        &coset_factor,
                    )
                    .expect("some lde");
            });
        });

        let inner_size = 1 << (log_size / 2);
        let new_domain = Domain::<Fr>::new_for_size(inner_size).expect("some domain");
        let inner_omega = new_domain.generator.pow(&[inner_size]);
        let precomputed_twiddles = hodor::fft::strided_fft::utils::precompute_twiddle_factors(
            &inner_omega,
            inner_size as usize,
        );

        let bench_id = BenchmarkId::new(format!("lde-with-square-root-fft-{}", size), num_cpus);
        group.bench_function(bench_id, |b| {
            b.iter(|| {
                poly2
                    .clone()
                    .lde_using_square_root_fft(&worker, lde_factor, &precomputed_twiddles)
            });
        });
    }
}
