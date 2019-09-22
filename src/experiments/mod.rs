pub(crate) mod vdf;
pub(crate) mod cubic_vdf;
mod square_root_calculator;

/*

This module contains a set of experiments like proving of a proper evaluation of VDF that itself is a square/cube root calculation.
Such case is ideal for Stark as it contains only trivial low-degree constraints that are dense and a very small number of registers

Interested user can find an example how to define the constraints using the convenience function in the `vdf.rs` or `cubic_vdf.rs` files

*/

mod tensor_lde;

use ff::*;

#[derive(PrimeField)]
#[PrimeFieldModulus = "3618502788666131213697322783095070105623107215331596699973092056135872020481"]
#[PrimeFieldGenerator = "3"]
pub struct Fr(FrRepr);

#[test]
fn test_proper_field(){
    let modulus = self::Fr::char();
    println!("Char = {}", modulus);
    let mul_gen = self::Fr::multiplicative_generator();
    let omega = Fr::root_of_unity();
    let one = self::Fr::one();
    let mut pow = self::Fr::one();
    println!("Max multiplicative subgroup of size 2^{}", Fr::S);
    for i in 0..Fr::S {
        let r = mul_gen.pow(&pow.into_repr());
        if r == one {
            panic!("R^(2^{}) = {}, but should not be one", i, r);
        }
        let o = omega.pow(&pow.into_repr());
        if o == one {
            panic!("Omega^(2^{}) = {}, but should not be one", i, r);
        }
        pow.double();
    }

    let o = omega.pow(&pow.into_repr());
    if o != one {
        panic!("R^(2^S) = {}, but should be one", o);
    }

    

}