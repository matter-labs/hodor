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

impl<F: PrimeField> std::fmt::Display for Constraint<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        writeln!(f, "Polynomial constraint of degree {}", self.degree)?;
        write!(f, "0 = {} ", self.constant_term)?;
        let num_terms = self.terms.len();
        for i in 0..num_terms {
            let term = &self.terms[i];
            write!(f, "{}", term)?;
            if i != num_terms - 1 {
                write!(f, " ")?;
            }
        }
        write!(f, ";")
    }
}

impl<F: PrimeField> std::fmt::Display for ConstraintTerm<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            ConstraintTerm::Univariate(t) => {
                write!(f, "{}", t)?;
            },
            ConstraintTerm::Polyvariate(p) => {
                let one = F::one();
                let mut minus_one = F::one();
                minus_one.negate();
                if p.coeff == one {
                    write!(f, "+ ")?;
                } else if p.coeff == minus_one {
                    write!(f, "- ")?;
                } else {
                    write!(f, "+ {}* ", p.coeff)?;
                }


                let num_terms = p.terms.len();
                for i in 0..num_terms {
                    let t = p.terms[i];
                    let num_steps = match t.steps_difference {
                        StepDifference::Steps(num_steps) => num_steps,
                        _ => unimplemented!()
                    };
                    match t.register {
                        Register::Register(reg_num) => {
                            write!(f, "(R_{}(t+{}))^{}", reg_num, num_steps, t.power)?;
                        },
                        _ => {
                            unimplemented!();
                        }
                    }
                    if i != num_terms - 1 {
                        write!(f, "*")?;
                    }
                }
            }
        }

        Ok(())
    }
}

impl<F: PrimeField> std::fmt::Display for UnivariateTerm<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        let one = F::one();
        let mut minus_one = F::one();
        minus_one.negate();
        if self.coeff == one {
            write!(f, "+ ")?;
        } else if self.coeff == minus_one {
            write!(f, "- ")?;
        } else {
            write!(f, "+ {}*", self.coeff)?;
        }
        let num_steps = match self.steps_difference {
            StepDifference::Steps(num_steps) => num_steps,
            _ => unimplemented!()
        };
        match self.register {
            Register::Register(reg_num) => {
                write!(f, "(R_{}(t+{}))^{}", reg_num, num_steps, self.power)?;
            },
            _ => {
                unimplemented!();
            }
        }

        Ok(())
    }
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

// Univariate term implementation

impl<F: PrimeField> From<Register> for UnivariateTerm<F> {
    fn from(register: Register) -> UnivariateTerm<F> {
        UnivariateTerm::<F> {
            coeff: F::one(),
            register: register,
            steps_difference: StepDifference::Steps(0),
            power: 1
        }
    }
}

impl<F: PrimeField> MulAssign<&F> for UnivariateTerm<F> {
    fn mul_assign(&mut self, rhs: &F) {
        self.coeff.mul_assign(&rhs);
    }
}

impl<F: PrimeField> MulAssign<F> for UnivariateTerm<F> {
    fn mul_assign(&mut self, rhs: F) {
        self.coeff.mul_assign(&rhs);
    }
}

impl<F: PrimeField> UnivariateTerm<F> {
    pub fn set_step_difference(&mut self, steps: usize) {
        self.steps_difference = StepDifference::Steps(steps);
    }

    pub fn pow(&self, power: u64) -> Self {
        let new_power = self.power.checked_mul(power).expect("maximum power is limited to U64MAX");
        let mut new = self.clone();
        new.power = new_power;

        new
    }
}

// Polyvariate term implementation

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

impl<F: PrimeField> MulAssign<UnivariateTerm<F>> for PolyvariateTerm<F> {
    fn mul_assign(&mut self, rhs: UnivariateTerm<F>) {
        let mut other = rhs;
        self.coeff.mul_assign(&other.coeff);
        other.coeff = F::one();
        self.total_degree += other.power;
        self.terms.push(other);
    }
}

impl<F: PrimeField> MulAssign<F> for PolyvariateTerm<F> {
    fn mul_assign(&mut self, other: F) {
        self.coeff.mul_assign(&other);
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

// Constraint implementation

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