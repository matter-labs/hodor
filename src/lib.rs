extern crate ff;
extern crate byteorder;
extern crate rand;
extern crate hex;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate blake2_rfc;
extern crate tiny_keccak;
extern crate crypto;

mod air;
pub mod fri;
pub mod utils;
pub mod fft;

use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};

#[derive(PrimeField)]
#[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
#[PrimeFieldGenerator = "7"]
pub struct Fr(FrRepr);