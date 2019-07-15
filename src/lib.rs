extern crate ff;
extern crate byteorder;
extern crate rand;
extern crate hex;
extern crate serde;
extern crate serde_derive;
extern crate blake2_rfc;
extern crate tiny_keccak;
extern crate crypto;

pub mod air;
pub mod arp;
pub mod fri;
pub mod utils;
pub mod fft;
pub mod domains;
pub mod polynomials;
pub mod ali;

use ff::{Field, PrimeField, PrimeFieldRepr};

// #[derive(PrimeField)]
// #[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
// #[PrimeFieldGenerator = "7"]
// pub struct Fr(FrRepr);

#[derive(PrimeField)]
#[PrimeFieldModulus = "257"]
#[PrimeFieldGenerator = "3"]
pub struct Fr(FrRepr);

#[derive(Debug)]
pub enum SynthesisError {
    Error,
}

use std::fmt;
use std::error::Error;

impl Error for SynthesisError {
    fn description(&self) -> &str {
        match *self {
            SynthesisError::Error => "General error for now",
        }
    }
}

impl fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "{}", self.description())
    }
}