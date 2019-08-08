use ff::{
    PrimeField,
};

use super::*;

use std::ops::{AddAssign, SubAssign, MulAssign};
use std::default::Default;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryConstraint<F: PrimeField> {
    pub register: Register,
    pub at_row: usize,
    pub value: Option<F>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Constraint<F: PrimeField> {
    pub constant_term: F,
    pub terms: Vec<ConstraintTerm<F>>,
    pub degree: u64,
    pub density: ConstraintDensity,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstraintTerm<F: PrimeField> {
    Univariate(UnivariateTerm<F>),
    Polyvariate(PolyvariateTerm<F>)
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct UnivariateTerm<F: PrimeField>{
    pub coeff: F,
    pub register: Register,
    pub steps_difference: StepDifference<F>,
    pub power: u64
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepDifference<F: PrimeField> {
    Steps(usize),
    Mask(F),
}

impl<F: PrimeField> std::hash::Hash for StepDifference<F> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match &self {
            StepDifference::Steps(s) => {
                s.hash(state);
            },
            StepDifference::Mask(m) => {
                m.into_raw_repr().as_ref().hash(state);
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyvariateTerm<F: PrimeField> {
    pub coeff: F,
    pub terms: Vec<UnivariateTerm<F>>,
    pub total_degree: u64
}

impl<F: PrimeField> From<(F, UnivariateTerm<F>)> for PolyvariateTerm<F> {
    fn from(w: (F, UnivariateTerm<F>)) -> PolyvariateTerm<F> {
        let (mut coeff, mut constraint) = w;
        // optimize
        coeff.mul_assign(&constraint.coeff);
        constraint.coeff = F::one();

        let degree = constraint.power;

        PolyvariateTerm::<F> {
            coeff: coeff,
            terms: vec![constraint],
            total_degree: degree
        }
    }
}

impl<F: PrimeField> Default for PolyvariateTerm<F> {
    fn default() -> PolyvariateTerm<F> {
        PolyvariateTerm::<F> {
            coeff: F::one(),
            terms: vec![],
            total_degree: 0u64
        }
    }
}

impl<F: PrimeField> Default for Constraint<F> {
    fn default() -> Constraint<F> {
        Constraint::<F> {
            constant_term: F::zero(),
            terms: vec![],
            degree: 0u64,
            density: ConstraintDensity::Dense(DenseConstraint::default()),
        }
    }
}

impl<F: PrimeField> MulAssign<UnivariateTerm<F>> for PolyvariateTerm<F> {
    fn mul_assign(&mut self, rhs: UnivariateTerm<F>) {
        let mut other = rhs;
        self.coeff.mul_assign(&other.coeff);
        other.coeff = F::one();
        self.total_degree += other.power;
        self.terms.push(other);
    }
}

impl<F: PrimeField> MulAssign<&F> for PolyvariateTerm<F> {
    fn mul_assign(&mut self, other: &F) {
        self.coeff.mul_assign(other);
    }
}

impl<F: PrimeField> MulAssign<PolyvariateTerm<F>> for PolyvariateTerm<F> {
    fn mul_assign(&mut self, rhs: PolyvariateTerm<F>) {
        self.coeff.mul_assign(&rhs.coeff);
        self.total_degree += rhs.total_degree;
        self.terms.extend(rhs.terms.into_iter());
    }
}

impl<F:PrimeField> AddAssign<PolyvariateTerm<F>> for Constraint<F> {
    fn add_assign(&mut self, rhs: PolyvariateTerm<F>) {
        if self.degree < rhs.total_degree {
            self.degree = rhs.total_degree;
        }
        self.terms.push(ConstraintTerm::Polyvariate(rhs));
    }
}

impl<F:PrimeField> SubAssign<PolyvariateTerm<F>> for Constraint<F> {
    fn sub_assign(&mut self, rhs: PolyvariateTerm<F>) {
        let mut rhs = rhs;
        rhs.coeff.negate();
        if self.degree < rhs.total_degree {
            self.degree = rhs.total_degree;
        }
        self.terms.push(ConstraintTerm::Polyvariate(rhs));
    }
}

impl<F:PrimeField> AddAssign<UnivariateTerm<F>> for Constraint<F> {
    fn add_assign(&mut self, rhs: UnivariateTerm<F>) {
        if self.degree < rhs.power {
            self.degree = rhs.power;
        }
        self.terms.push(ConstraintTerm::Univariate(rhs));
    }
}

impl<F:PrimeField> SubAssign<UnivariateTerm<F>> for Constraint<F> {
    fn sub_assign(&mut self, rhs: UnivariateTerm<F>) {
        if self.degree < rhs.power {
            self.degree = rhs.power;
        }
        let mut other = rhs;
        other.coeff.negate();
        self.terms.push(ConstraintTerm::Univariate(other));
    }
}

impl<F:PrimeField> AddAssign<&F> for Constraint<F> {
    fn add_assign(&mut self, rhs: &F) {
        self.constant_term.add_assign(rhs);
    }
}

impl<F:PrimeField> SubAssign<&F> for Constraint<F> {
    fn sub_assign(&mut self, rhs: &F) {
        self.constant_term.sub_assign(rhs);
    }
}