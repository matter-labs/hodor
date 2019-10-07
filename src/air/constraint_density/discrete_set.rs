use super::*;
use crate::ff::PrimeField;
// Happens at the subdomain (e.g. every 4th row) of the column. May be displsced
// to happen on e.g. 4k + 1
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct DiscreteSetConstraint {
    pub rows: Vec<usize>
}

pub(crate) struct DiscreteSetConstraintQuery {
    rows: Vec<usize>,
    current: usize,
}

impl DensityQuery for DiscreteSetConstraintQuery {
    fn next_row(&mut self) -> Option<usize> {
        if self.current < self.rows.len() {
            let t = self.rows[self.current];

            self.current += 1;

            Some(t)
        } else {
            None
        }
    }
}

impl DiscreteSetConstraint {
    pub fn new(mut rows: Vec<usize>) -> Self {
        rows.sort();

        Self {
            rows
        }
    }

    fn roots<F: PrimeField>(&self, column_domain: &Domain<F>) -> Vec<F> {
        let roots_generator = column_domain.generator;

        let mut rows = self.rows.clone();
        rows.sort();

        let mut roots = vec![];
        for num in rows.iter() {
            let x = roots_generator.pow([*num as u64]);
            roots.push(x);               
        }

        roots
    }
}

impl<F: PrimeField> ConstraintDensity<F> for DiscreteSetConstraint {
    fn calculate_density_query(&self, num_rows_in_trace: usize) -> Box<dyn DensityQuery> {
        assert!(self.rows.len() > 0);
        let max = *self.rows.iter().max().unwrap();
        assert!(max < num_rows_in_trace);

        let mut rows = self.rows.clone();
        rows.sort();

        let query = DiscreteSetConstraintQuery {
            rows: rows,
            current: 0,
        };

        Box::from(query)
    }

    fn inverse_divisor_in_coset(   
        &self,         
        column_domain: &Domain<F>,
        evaluation_domain: &Domain<F>,
        precomputations_for_column_domain: &Option<PrecomputedOmegas<F>>,
        precomputations_for_evaluation_domain: &Option<PrecomputedOmegas<F>>,
        num_rows_in_trace: usize,
        worker: &Worker
    ) -> Result<(Polynomial<F, Values>, usize), SynthesisError> {
        let roots = self.roots(&column_domain);

        let mut poly = Polynomial::<F, Coefficients>::from_roots(roots, &worker)?;

        let divisor_degree = self.rows.len();

        if poly.size() < worker.cpus {
            poly.pad_to_size(worker.cpus)?;
        }

        let factor = (evaluation_domain.size as usize) / poly.size();

        let mut inv_divisors = poly.coset_lde(&worker, factor)?;

        inv_divisors.batch_inversion(&worker)?;

        Ok((inv_divisors, divisor_degree))
    }

    fn evaluate_inversed_at(
        &self,            
        at: F,
        column_domain: &Domain<F>,
        evaluation_domain: &Domain<F>,
        precomputations_for_column_domain: &Option<PrecomputedOmegas<F>>,
        precomputations_for_evaluation_domain: &Option<PrecomputedOmegas<F>>,
        num_rows_in_trace: usize,
    ) -> Result<F, SynthesisError> {
        let roots = self.roots(&column_domain);

        let mut result = F::one();

        for root in roots.into_iter() {
            let mut tmp = at;
            tmp.sub_assign(&root);
            result.mul_assign(&tmp);
        }

        let inv_vanishing_at_subset = result.inverse().unwrap();

        Ok(inv_vanishing_at_subset)
    }

    fn is_equal_to(&self, other: &(dyn Any + Send + Sync)) -> bool {
        if let Some(other) = other.downcast_ref::<Self>() {
            self == other
        } else {
            false
        }
    }

    fn as_any(&self) -> &(dyn Any + Send + Sync) {
        self
    }
}