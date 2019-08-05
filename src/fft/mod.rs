pub(crate) mod multicore;
pub(crate) mod fft;
pub(crate) mod lde;

// pub(crate) mod recursive_fft;
// pub(crate) mod recursive_lde;

use cfg_if;

#[cfg(feature = "nightly")]
mod prefetch_lde;
#[cfg(feature = "nightly")]
mod prefetch_fft;
#[cfg(feature = "nightly")]
mod prefetch;

use ff::PrimeField;
use self::multicore::Worker;

cfg_if! {
    if #[cfg(feature = "nightly")] {
        #[inline(always)]
        pub(crate) fn best_lde<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, lde_factor: usize) {
            self::prefetch_lde::best_lde(a, worker, omega, log_n, lde_factor)
        }

        #[inline(always)]
        pub(crate) fn best_fft<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, use_cpus_hint: Option<usize>) {
            self::prefetch_fft::best_fft(a, worker, omega, log_n, use_cpus_hint)
        }

        #[inline(always)]
        pub(crate) fn serial_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32) {
            self::prefetch_fft::serial_fft(a, omega, log_n)
        }
    } else {
        #[inline(always)]
        pub(crate) fn best_lde<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, lde_factor: usize) {
            self::lde::best_lde(a, worker, omega, log_n, lde_factor)
        }
        #[inline(always)]
        pub(crate) fn best_fft<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, use_cpus_hint: Option<usize>) {
            self::fft::best_fft(a, worker, omega, log_n, use_cpus_hint)
        }
        #[inline(always)]
        pub(crate) fn serial_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32) {
            self::fft::serial_fft(a, omega, log_n)
        }
    }  
}


pub(crate) mod dit_fft;
pub(crate) mod radix4_fft;


#[test]
fn test_sequential_FFT()
{
    use rand::{XorShiftRng, SeedableRng, Rand};
    const LOG_N: u32 = 10;
    const N: usize = 1 << LOG_N;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    use ff::Field;
    use crate::experiments::vdf::Fr;
    use std::time::Instant;
    use crate::domains::Domain;

    let mut a = (0..N).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let mut b = a.clone();
    let mut c = a.clone();

    let domain = Domain::<Fr>::new_for_size(a.len() as u64).unwrap();
    let omega = domain.generator;

    let mut start = Instant::now();
    fft::serial_fft::<Fr>(&mut a, &omega, LOG_N);
    let mut end = Instant::now();
    let radix_2_time = end - start;

    let mut start = Instant::now();
    radix4_fft::serial_fft_radix_4::<Fr>(&mut b, &omega, LOG_N);
    let mut end = Instant::now();
    let radix_4_time = end - start;

    start = Instant::now();
    dit_fft::serial_DIT_fft::<Fr>(&mut c, &omega, LOG_N, N);
    end = Instant::now();
    let dit_fft_time = end - start;

    println!("Radix2 time: {}", radix_2_time.subsec_millis());
    println!("Radix4 time: {}", radix_4_time.subsec_millis());
    println!("dit time: {}", dit_fft_time.subsec_millis());
    

    // println!("{:?}", a);
    // println!("{:?}", b);

    let matching_radix4 = a.iter().zip(b.iter()).filter(|(a, b)| *a == *b).count();
    let matching_dit = a.iter().zip(c.iter()).filter(|(a, c)| *a == *c).count();
    assert_eq!(matching_radix4, N);
    assert_eq!(matching_dit, N);
    
}

pub fn distribute_powers<F: PrimeField>(coeffs: &mut [F], worker: &Worker, g: F)
{
    worker.scope(coeffs.len(), |scope, chunk| {
        for (i, v) in coeffs.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut u = g.pow(&[(i * chunk) as u64]);
                for v in v.iter_mut() {
                    v.mul_assign(&u);
                    u.mul_assign(&g);
                }
            });
        }
    });
}



#[test]
fn test_parallel_FFT()
{
    use rand::{XorShiftRng, SeedableRng, Rand};
    const LOG_N: u32 = 10;
    const N: usize = 1 << LOG_N;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    use ff::Field;
    use crate::experiments::vdf::Fr;
    use std::time::Instant;
    use crate::domains::Domain;
    use crate::fft::multicore::Worker;

    //create two different workers: general and special worker for radix4 FFT
    let general_worker = Worker::new();
    let log_cpus = general_worker.log_num_cpus();
    let radix4_log_cpus = if log_cpus % 2 == 0 {log_cpus} else {log_cpus - 1};
    let radix4_worker = Worker::new_with_cpus(1 << radix4_log_cpus);

    let mut a = (0..N).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let mut b = a.clone();
    let mut c = a.clone();

    let domain = Domain::<Fr>::new_for_size(a.len() as u64).unwrap();
    let omega = domain.generator;

    let mut start = Instant::now();
    fft::parallel_fft::<Fr>(&mut a, &general_worker, &omega, LOG_N, log_cpus);
    let mut end = Instant::now();
    let radix_2_time = end - start;

    let mut start = Instant::now();
    radix4_fft::parallel_fft_radix_4::<Fr>(&mut b, &radix4_worker, &omega, LOG_N, radix4_log_cpus);
    let mut end = Instant::now();
    let radix_4_time = end - start;

    start = Instant::now();
    dit_fft::parallel_DIT_fft::<Fr>(&mut c, &general_worker, &omega, LOG_N, log_cpus, N);
    end = Instant::now();
    let dit_fft_time = end - start;

    println!("Radix2 time: {}", radix_2_time.subsec_millis());
    println!("Radix4 time: {}", radix_4_time.subsec_millis());
    println!("dit time: {}", dit_fft_time.subsec_millis());
    
    let matching_radix4 = a.iter().zip(b.iter()).filter(|(a, b)| *a == *b).count();
    let matching_dit = a.iter().zip(c.iter()).filter(|(x, y)| *x == *y).count();

    // println!("{:?}", a);
    // println!("{:?}", b);
    
    assert_eq!(matching_radix4, N);
    assert_eq!(matching_dit, N);
    
}


#[test]
fn test_FFT_Prunning()
{
    use rand::{XorShiftRng, SeedableRng, Rand};
    const LOG_NONZERO_N: u32 = 18;
    const LOG_N: u32 = LOG_NONZERO_N + 4;
    const N: usize = 1 << LOG_N;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    use ff::Field;
    use crate::experiments::vdf::Fr;
    use std::time::Instant;
    use crate::domains::Domain;
    use crate::fft::multicore::Worker;

    let worker = Worker::new();

    let mut a = (0..N).map(|i| if i < (1<<LOG_NONZERO_N) {Fr::rand(rng)} else {Fr::zero()}).collect::<Vec<_>>();
    let mut b = a.clone();

    let domain = Domain::<Fr>::new_for_size(a.len() as u64).unwrap();
    let omega = domain.generator;

    let mut start = Instant::now();
    //fft::parallel_fft::<Fr>(&mut a, &worker, &omega, LOG_N, worker.log_num_cpus());
    dit_fft::parallel_DIT_fft::<Fr>(&mut a, &worker, &omega, LOG_N, worker.log_num_cpus(), N);
    let mut end = Instant::now();
    let usual_time = end - start;

    start = Instant::now();
    dit_fft::parallel_DIT_fft::<Fr>(&mut b, &worker, &omega, LOG_N, worker.log_num_cpus(), 1<<LOG_NONZERO_N);
    end = Instant::now();
    let pruning_time = end - start;

    println!("usual time: {}", usual_time.subsec_millis());
    println!("pruning time: {}", pruning_time.subsec_millis());
    
    let matching = a.iter().zip(b.iter()).filter(|(a, b)| *a == *b).count();

    // println!("{:?}", a);
    // println!("{:?}", b);
    
    assert_eq!(matching, N); 
}


// #[cfg(feature = "nightly")]
// extern crate prefetch;

// #[cfg(feature = "nightly")]
// pub(crate) mod radix4_fft;

// #[test]
// fn test_large_lde() {
//     use rand::{XorShiftRng, SeedableRng, Rand};
//     const LOG_N: usize = 22;
//     const BASE: usize = 1 << LOG_N;
//     const LOG_LDE: usize = 7;
//     const LDE_FACTOR: usize = 1 << LOG_LDE;
//     let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

//     use ff::Field;
//     use crate::experiments::vdf::Fr;
//     use crate::domains::Domain;
//     use crate::fft::multicore::Worker;
//     use crate::polynomials::Polynomial;
//     use std::time::Instant;

//     let worker = Worker::new();

//     let mut coeffs = (0..BASE).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

//     let mut poly = Polynomial::from_coeffs(coeffs.clone()).unwrap();

//     poly.pad_by_factor(LDE_FACTOR).unwrap();
//     let now = Instant::now();
//     let naive_lde = poly.fft(&worker);
//     println!("naive LDE taken {}ms", now.elapsed().as_millis());

//     coeffs.resize(BASE * LDE_FACTOR, Fr::zero());

//     let domain = Domain::<Fr>::new_for_size(coeffs.len() as u64).unwrap();
//     let omega = domain.generator;
//     let log_n = (LOG_N + LOG_LDE) as u32;

//     let now = Instant::now();
//     let mut lde = coeffs.clone();
//     self::best_lde(&mut lde, &worker, &omega, log_n, LDE_FACTOR);
//     println!("LDE taken {}ms", now.elapsed().as_millis());

//     assert!(naive_lde.into_coeffs() == lde);
// }