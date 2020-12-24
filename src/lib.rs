#![feature(min_const_generics, asm)]
#![allow(dead_code)]


extern crate ff;
extern crate byteorder;
extern crate rand;
extern crate hex;
extern crate tiny_keccak;
extern crate indexmap;
extern crate crypto;
extern crate blake2b_simd;
extern crate blake2s_simd;
#[macro_use] extern crate lazy_static;
#[macro_use] extern crate cfg_if;

pub mod air;
pub mod arp;
pub mod fri;
pub mod utils;
pub mod fft;
pub mod domains;
pub mod polynomials;
pub mod ali;
pub mod iop;
pub mod transcript;
pub mod precomputations;
pub mod verifier;
pub mod prover;

pub mod optimized_fields;

mod experiments;

pub(crate) mod bn256;
pub mod f125;

use ff::{Field, PrimeField, PrimeFieldRepr};

#[derive(PrimeField)]
#[PrimeFieldModulus = "257"]
#[PrimeFieldGenerator = "3"]
pub struct Fr(FrRepr);

#[derive(Debug)]
pub enum SynthesisError {
    Error,
    Unsatisfied(String),
    InvalidValue(String),
    DivisionByZero(String),
}

use std::fmt;
use std::error::Error;

impl Error for SynthesisError {}

impl fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            SynthesisError::Error => write!(f, "{}", "General error for now"),
            SynthesisError::Unsatisfied(descr) => write!(f, "Unsatisfied constraint, {}", descr),
            SynthesisError::InvalidValue(descr) => write!(f, "Invalid parameter value, {}", descr),
            SynthesisError::DivisionByZero(descr) => write!(f, "Division by zero, {}", descr),
        }
    }
}