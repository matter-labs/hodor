extern crate ff;
extern crate byteorder;
extern crate rand;
extern crate hex;
extern crate serde;
#[macro_use]
extern crate serde_derive;

mod field;
mod air;
mod fri;

use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};

#[derive(PrimeField)]
#[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
#[PrimeFieldGenerator = "7"]
pub struct Fr(FrRepr);