#![feature(asm, test)]

use criterion::Criterion;
mod transposition;
mod arith;
mod fft;
mod comparisons;

fn main() {
    let crit = &mut Criterion::default().configure_from_args();
    // arith::group(crit);
    // transposition::group(crit);
    // fft::group(crit);
    comparisons::group(crit);
    crit.final_summary();
}