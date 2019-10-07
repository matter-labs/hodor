use super::*;

// Happens at the subdomain (e.g. every 4th row) of the column. May be displsced
// to happen on e.g. 4k + 1
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct SubdomainConstraint {
    pub subdomain_factor: usize,
    pub subdomain_coset_number: usize,
}

pub(crate) struct SubdomainConstraintQuery {
    subdomain_factor: usize,
    subdomain_coset_number: usize,
    subdomain_element: usize,
    subdomain_size: usize,
    limit: usize
}

impl DensityQuery for SubdomainConstraintQuery {
    fn next_row(&mut self) -> Option<usize> {
        if self.subdomain_element < self.subdomain_size {
            let t = self.subdomain_element;
            let t = t * self.subdomain_factor + self.subdomain_coset_number;

            self.subdomain_element += 1;

            Some(t)
        } else {
            None
        }
    }
}

impl<F: PrimeField> ConstraintDensity<F> for SubdomainConstraint {
    fn calculate_density_query(&self, num_rows_in_trace: usize) -> Box<dyn DensityQuery> {
        assert!(num_rows_in_trace % self.subdomain_factor == 0);
        assert!(self.subdomain_coset_number < self.subdomain_factor);

        let query = SubdomainConstraintQuery {
            subdomain_factor: self.subdomain_factor,
            subdomain_coset_number: self.subdomain_coset_number,
            subdomain_element: 0,
            subdomain_size: num_rows_in_trace / self.subdomain_factor,
            limit: 0
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
        assert!(self.subdomain_coset_number < self.subdomain_factor);
        assert!((column_domain.size as usize) % self.subdomain_factor == 0);

        let divisor_degree = (column_domain.size as usize) / self.subdomain_factor;

        let column_domain_generator = column_domain.generator;
        let evaluation_domain_generator = evaluation_domain.generator;
        let multiplicative_generator = F::multiplicative_generator();

        let subset_size = divisor_degree;

        let evaluation_domain_size = evaluation_domain.size as usize;
        // assert!(evaluation_domain_size == precomputations.coset.len());

        // updated strategy:
        // - calculate {X}^T - 1 in the coset of the column domain. This is a single powering and inversion
        // - multiply by (X - root) in the coset for all the roots
        // - coset IFFT and then coset LDE

        let subset_domain = Domain::<F>::new_for_size(subset_size as u64)?;

        let mut vanishing_at_subset = F::multiplicative_generator().pow([subset_size as u64]);
        let coset_factor = subset_domain.generator.pow([self.subdomain_coset_number as u64]);
        vanishing_at_subset.sub_assign(&coset_factor);

        let inv_vanishing_at_subset = vanishing_at_subset.inverse().unwrap();

        let coset_values = vec![inv_vanishing_at_subset; subset_size];

        let inv_divisors = Polynomial::<F, Values>::from_values(coset_values)?;
        let inv_divisors = inv_divisors.icoset_fft(&worker);

        let factor = (evaluation_domain.size as usize) / subset_size;

        let mut inv_divisors = inv_divisors.coset_lde(&worker, factor)?;

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

        let subset_size = (column_domain.size as usize) / self.subdomain_factor;

        let subset_domain = Domain::<F>::new_for_size(subset_size as u64)?;
        
        let mut vanishing_at_subset = at.pow([subset_size as u64]);
        let coset_factor = subset_domain.generator.pow([self.subdomain_coset_number as u64]);
        vanishing_at_subset.sub_assign(&coset_factor);

        let inv_vanishing_at_subset = vanishing_at_subset.inverse().unwrap();

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