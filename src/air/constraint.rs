use ff::{
    PrimeField,
};

use super::Register;

use std::ops::{AddAssign, SubAssign, MulAssign};
use std::default::Default;

#[derive(Debug, Clone)]
pub struct UnivariateConstraintTerm<F: PrimeField>(pub F, pub (Register, usize), pub u64); // coeff, (register, step delta), power

#[derive(Debug, Clone)]
pub struct PolyvariateConstraintTerm<F: PrimeField>(pub F, pub Vec<UnivariateConstraintTerm<F>>); // coeff, terms

#[derive(Debug, Clone)]
pub struct PolynomialConstraint<F: PrimeField>(pub F, pub Vec<PolynomialConstraintTerm<F>>); // constant and sum of terms

#[derive(Debug, Clone)]
pub enum PolynomialConstraintTerm<F: PrimeField> {
    Univariate(UnivariateConstraintTerm<F>),
    Polyvariate(PolyvariateConstraintTerm<F>)
}

impl<F: PrimeField> From<(F, UnivariateConstraintTerm<F>)> for PolyvariateConstraintTerm<F> {
    fn from(w: (F, UnivariateConstraintTerm<F>)) -> PolyvariateConstraintTerm<F> {
        PolyvariateConstraintTerm(w.0, vec![w.1])
    }
}

impl<F: PrimeField> Default for PolyvariateConstraintTerm<F> {
    fn default() -> PolyvariateConstraintTerm<F> {
        PolyvariateConstraintTerm(F::one(), vec![])
    }
}

impl<F: PrimeField> Default for PolynomialConstraint<F> {
    fn default() -> PolynomialConstraint<F> {
        PolynomialConstraint(F::zero(), vec![])
    }
}

impl<F: PrimeField> MulAssign<UnivariateConstraintTerm<F>> for PolyvariateConstraintTerm<F> {
    fn mul_assign(&mut self, other: UnivariateConstraintTerm<F>) {
        self.1.push(other);
    }
}

impl<F: PrimeField> MulAssign<&F> for PolyvariateConstraintTerm<F> {
    fn mul_assign(&mut self, other: &F) {
        self.0.mul_assign(other);
    }
}

impl<F: PrimeField> MulAssign<PolyvariateConstraintTerm<F>> for PolyvariateConstraintTerm<F> {
    fn mul_assign(&mut self, rhs: PolyvariateConstraintTerm<F>) {
        self.0.mul_assign(&rhs.0);
        self.1.extend(rhs.1.into_iter());
    }
}

impl<F:PrimeField> AddAssign<PolyvariateConstraintTerm<F>> for PolynomialConstraint<F> {
    fn add_assign(&mut self, rhs: PolyvariateConstraintTerm<F>) {
        self.1.push(PolynomialConstraintTerm::Polyvariate(rhs));
    }
}

impl<F:PrimeField> AddAssign<UnivariateConstraintTerm<F>> for PolynomialConstraint<F> {
    fn add_assign(&mut self, rhs: UnivariateConstraintTerm<F>) {
        self.1.push(PolynomialConstraintTerm::Univariate(rhs));
    }
}

impl<F:PrimeField> SubAssign<UnivariateConstraintTerm<F>> for PolynomialConstraint<F> {
    fn sub_assign(&mut self, rhs: UnivariateConstraintTerm<F>) {
        let mut other = rhs;
        other.0.negate();
        self.1.push(PolynomialConstraintTerm::Univariate(other));
    }
}

impl<F:PrimeField> AddAssign<&F> for PolynomialConstraint<F> {
    fn add_assign(&mut self, rhs: &F) {
        self.0.add_assign(rhs);
    }
}

impl<F:PrimeField> SubAssign<&F> for PolynomialConstraint<F> {
    fn sub_assign(&mut self, rhs: &F) {
        self.0.sub_assign(rhs);
    }
}