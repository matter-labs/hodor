use ff::PrimeField;

use crate::air::*;
use crate::domains::*;

pub(crate) fn remap_univariate_term<F: PrimeField>(
    term: &mut UnivariateTerm<F>,
    column_domain: &Domain<F>,
) {
    let step_delta = match term.steps_difference {
        StepDifference::Steps(num) => {
            num as u64
        },
        _ => {
            unreachable!("Step differences are not masks yet");
        }
    };

    let mask = column_domain.generator.pow([step_delta]);

    term.steps_difference = StepDifference::Mask(mask);
}

pub(crate) fn remap_term<F: PrimeField>(
    term: &mut ConstraintTerm<F>,
    column_domain: &Domain<F>,
) {
    match term {
        ConstraintTerm::Univariate(ref mut t) => {
            remap_univariate_term(
                t, 
                &column_domain
            );
        },
        ConstraintTerm::Polyvariate(ref mut poly_term) => {
            for mut t in poly_term.terms.iter_mut() {
                remap_univariate_term(
                    &mut t, 
                    &column_domain
                );
            }
        }
    }
}

pub(crate) fn remap_constraint<F: PrimeField>(
    constraint: &mut Constraint<F>,
    column_domain: &Domain<F>,
) {
    for mut t in constraint.terms.iter_mut() {
        remap_term(
            &mut t, 
            &column_domain
        );
    }
}