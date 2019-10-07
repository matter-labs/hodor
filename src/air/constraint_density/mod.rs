mod density_query;

pub use self::density_query::DensityQuery;

use crate::ff::PrimeField;
use crate::domains::Domain;
use crate::precomputations::PrecomputedOmegas;
use crate::fft::multicore::Worker;
use crate::SynthesisError;
use crate::polynomials::*;
use std::any::Any;

mod dense_constraint;
mod subdomain_constraint;
mod discrete_set;
mod anywhere_but_in_set;

pub use dense_constraint::DenseConstraint;
pub use subdomain_constraint::SubdomainConstraint;
pub use discrete_set::DiscreteSetConstraint;
pub use anywhere_but_in_set::AnywhereButInSetConstraint;

pub trait ConstraintDensity<F: PrimeField>:
    objekt::Clone
    + Send 
    + Sync 
    + 'static 
    + Any
    + std::fmt::Debug
    // + std::hash::Hash
{
    fn calculate_density_query(&self, num_rows_in_trace: usize) -> Box<dyn DensityQuery>;

    fn inverse_divisor_in_coset(
        &self,            
        column_domain: &Domain<F>,
        evaluation_domain: &Domain<F>,
        precomputations_for_column_domain: &Option<PrecomputedOmegas<F>>,
        precomputations_for_evaluation_domain: &Option<PrecomputedOmegas<F>>,
        num_rows_in_trace: usize,
        worker: &Worker
    ) -> Result<(Polynomial<F, Values>, usize), SynthesisError>;

    fn evaluate_inversed_at(
        &self,            
        at: F,
        column_domain: &Domain<F>,
        evaluation_domain: &Domain<F>,
        precomputations_for_column_domain: &Option<PrecomputedOmegas<F>>,
        precomputations_for_evaluation_domain: &Option<PrecomputedOmegas<F>>,
        num_rows_in_trace: usize,
    ) -> Result<F, SynthesisError>;

    fn is_equal_to(&self, other: &(dyn Any + Send + Sync)) -> bool;

    fn as_any(&self) -> &(dyn Any + Send + Sync);
}

objekt::clone_trait_object!(<F: PrimeField> ConstraintDensity<F>);

impl<F: PrimeField> PartialEq<dyn ConstraintDensity<F>> for Box<dyn ConstraintDensity<F>> {
    fn eq(&self, other: &dyn ConstraintDensity<F>) -> bool {
        self.is_equal_to(other.as_any())
    }
}

impl<F: PrimeField> PartialEq for dyn ConstraintDensity<F> {
    fn eq(&self, other: &dyn ConstraintDensity<F>) -> bool {
        self.is_equal_to(other.as_any())
    }
}

impl<F: PrimeField> Eq for dyn ConstraintDensity<F> {}

impl<F: PrimeField> std::hash::Hash for Box<dyn ConstraintDensity<F>> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.as_any().type_id().hash(state);
    }
}