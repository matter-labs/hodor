use criterion::{AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use rand::{XorShiftRng, SeedableRng, Rand, Rng};
use hodor::{domains::Domain, optimized_fields::f252_asm_generated::Fr, polynomials::Polynomial};
use hodor::fft::multicore::Worker;
use hodor::ff::{PrimeField, Field};

pub fn group(crit: &mut Criterion){    
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    // let worker = Worker::new();

    let mut group = crit.benchmark_group("LDE");
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);


    for num_cpus in (2..50).step_by(2){
        let worker = Worker::new_with_cpus(num_cpus);
        let log_size = 22;
        let size = 1<<log_size;
        let lde_factor = 2;
        let coeffs: Vec<_> = (0..size).map(|_| Fr::rand(rng)).collect();    
        let poly1 = Polynomial::from_coeffs(coeffs).expect("some poly");
        let poly2 = poly1.clone();
    
        let bench_id = BenchmarkId::new(format!("lde-with-best-fft-{}", size), num_cpus);
        group.bench_function(bench_id, |b|{
            b.iter(||{
                poly1.clone().lde(&worker, lde_factor).expect("some lde");
            });
        });

        let inner_size = 1<< (log_size/2);
        let new_domain = Domain::<Fr>::new_for_size(inner_size).expect("some domain");
        let inner_omega = new_domain.generator.pow(&[inner_size]);
        let precomputed_twiddles = hodor::fft::strided_fft::utils::precompute_twiddle_factors(&inner_omega, inner_size as usize);
        
        let bench_id = BenchmarkId::new(format!("lde-with-square-root-fft-{}", size), num_cpus);
        group.bench_function(bench_id, |b|{
            b.iter(||{
                poly2.clone().lde_using_square_root_fft(&worker, lde_factor, &precomputed_twiddles)
            });
        });
    }
    
}