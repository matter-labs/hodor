use crate::ff::*;

#[derive(PrimeField)]
#[PrimeFieldModulus = "63802944035360449460622495747797942273"]
#[PrimeFieldGenerator = "3"]
pub struct Fr(FrRepr);


pub fn dummy(a: u64, b: u64) -> String {
    let mut a = Fr::from_str(&a.to_string()).unwrap();
    let b = Fr::from_str(&b.to_string()).unwrap();

    a.mul_assign(&b);

    format!("{}", a)
}