pub mod fp2;

use super::Fr;
use ff::*;

#[test]
fn calculate() {
    let three = Fr::from_str("3").unwrap();
    let three_repr = three.into_raw_repr();
    for e in three_repr.as_ref().iter(){
        println!("0x{:016x}", e);
    }
    
}