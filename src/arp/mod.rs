use crate::air::{Register, PolynomialConstraint, ConstraintDensity};
use ff::PrimeField;

pub trait IntoARP<F: PrimeField> {
    // return full trace, trace constraints and boundary constraints
    fn into_arp(self) -> (Vec<Vec<F>>, (PolynomialConstraint<F>, ConstraintDensity), (Register, usize, F));
}