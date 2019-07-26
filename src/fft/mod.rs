pub(crate) mod multicore;
pub(crate) mod fft;
pub mod radix2_domain;

pub(crate) mod recursive_fft;
pub(crate) mod recursive_lde;

use cfg_if;

#[cfg(feature = "nightly")]
pub(crate) mod prefetch_lde;

#[cfg(not(feature = "nightly"))]
pub(crate) mod lde;

use ff::PrimeField;
use self::multicore::Worker;

cfg_if! {
    if #[cfg(feature = "nightly")] {
        #[inline(always)]
        pub(crate) fn best_lde<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, lde_factor: usize) {
            self::prefetch_lde::best_lde(a, worker, omega, log_n, lde_factor)
        }
    } else {
        #[inline(always)]
        pub(crate) fn best_lde<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, lde_factor: usize) {
            self::lde::best_lde(a, worker, omega, log_n, lde_factor)
        }
    }  
}

pub(crate) mod dit_fft;
pub(crate) mod radix4_fft;


#[test]
fn test_FFT()
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
    let matching_dit = a.iter().zip(b.iter()).filter(|(a, c)| *a == *c).count();
    assert_eq!(matching_radix4, N);
    assert_eq!(matching_dit, N);
    
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