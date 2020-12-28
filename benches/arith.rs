use std::convert::*;
use criterion::{black_box, Criterion, Bencher};
use rand::{Rng, XorShiftRng, SeedableRng};
use hodor::ff::*;
use criterion_convenience::*;

fn add_assing_256(c: &mut Criterion) {
    use hodor::optimized_fields::f252::Fr;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Add assign 256", |bencher| bencher.iter(|| black_box(a).add_assign(&black_box(b))));
}

fn add_assing_256_asm(c: &mut Criterion) {
    use hodor::optimized_fields::f252::Fr;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Add assign 256 ASM", |bencher| bencher.iter(|| black_box(a).add_assign(&black_box(b))));
}

fn mul_assing_128(c: &mut Criterion) {
    use hodor::optimized_fields::naive_f125::Fr;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Mont mul assign 128", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

fn mul_assing_128_asm(c: &mut Criterion) {
    use hodor::optimized_fields::f125::Fr;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let mut a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Mont mul assign 128 ASM manual", |bencher| bencher.iter(|| a.mul_assign(&b)));
    // c.bench_function("Mont mul assign 256 ASM", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

// fn mul_assing_256(c: &mut Criterion) {
//     use rand::{Rng, XorShiftRng, SeedableRng};
//     use fr256::Fr;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fr = rng.gen();
//     let b: Fr = rng.gen();

//     c.bench_function("Mont mul assign 256", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
// }

// fn mul_assing_256_asm_derive(c: &mut Criterion) {
//     use rand::{Rng, XorShiftRng, SeedableRng};
//     use fr256::Fr;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let mut a: Fr = rng.gen();
//     let b: Fr = rng.gen();

//     c.bench_function("Mont mul assign 256 ASM derived", |bencher| bencher.iter(|| a.mul_assign(&b)));
//     // c.bench_function("Mont mul assign 256", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
// }

// fn mul_assing_256_asm(c: &mut Criterion) {
//     use rand::{Rng, XorShiftRng, SeedableRng};
//     use hodor::optimized_fields::f252::Fr;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let mut a: Fr = rng.gen();
//     let b: Fr = rng.gen();

//     c.bench_function("Mont mul assign 256 ASM manual", |bencher| bencher.iter(|| a.mul_assign(&b)));
//     // c.bench_function("Mont mul assign 256 ASM", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
// }

// fn fft_rec_small(c: &mut Criterion) {
//     use rand::{Rng, XorShiftRng, SeedableRng};
//     use fr::Fr;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

//     const SIZE: usize = 1 << 12;
//     let domain = hodor::domains::Domain::<Fr>::new_for_size(SIZE as u64).unwrap();
//     let twiddles = hodor::fft::strided_fft::utils::precompute_twiddle_factors(&domain.generator, SIZE);
//     let values: Vec<Fr> = (0..SIZE).map(|_| rng.gen::<Fr>()).collect();
//     let mut values: [Fr; SIZE] = values.try_into().unwrap();
//     c.bench_function("Small recursive fft", 
//         |bencher| bencher.iter(
//             || hodor::fft::strided_fft::fft::small_size_serial_fft::<Fr, SIZE, 128>(&mut values, &twiddles, 0, 1, 1)
//         )
//     );
// }

pub fn group(crit: &mut Criterion) {
    mul_assing_128(crit);
    mul_assing_128_asm(crit);
}
