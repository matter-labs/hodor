pub mod f252;
pub mod naive_f252;
pub mod f125;
pub mod naive_f125;

use ff::*;

pub trait PrimeFieldExt: PrimeField {
    fn add(a: Self, b: Self) -> Self;
    fn sub(a: Self, b: Self) -> Self;
    fn mul(a: Self, b: Self) -> Self;
    fn minus_one() -> Self;
}

pub trait OneBitPartialRed: PrimeFieldExt {
    fn add_assign_noreduce(&mut self, b: &Self);
    fn sub_assign_and_add_modulus(&mut self, b: &Self);
    fn mul_assign_noreduce(&mut self, b: &Self);
    fn sub_modulus(&mut self);
    fn reduce_completely(&mut self);
}

pub trait TwoBitsPartialRed: OneBitPartialRed {
    fn sub_two_moduluses(&mut self);
}

pub trait ThreeBitsPartialRed: TwoBitsPartialRed {
    fn sub_four_moduluses(&mut self);
}


