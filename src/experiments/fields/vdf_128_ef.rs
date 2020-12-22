use super::*;

#[derive(PrimeField)]
#[PrimeFieldModulus = "63802944035360449460622495747797942273"]
#[PrimeFieldGenerator = "3"]
pub struct Fr(FrRepr);