use super::*;

// Happens at every row including start_at
// Span is an auxilary information that will not allow to "wrap around" the trace
// and also can be used to stop earlier
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct AnywhereButInSetConstraint {
    pub not_in_rows: Vec<usize>
}


impl AnywhereButInSetConstraint {
    pub fn new(mut not_in: Vec<usize>) -> Self {
        
        not_in.sort();

        AnywhereButInSetConstraint{
            not_in_rows: not_in,
        }
    }
}

impl std::default::Default for AnywhereButInSetConstraint {
    fn default() -> Self {
        AnywhereButInSetConstraint{
            not_in_rows: vec![]
        }
    }
}


pub(crate) struct AnywhereButInSetConstraintQuery {
    not_in_rows: Vec<usize>,
    current_row: usize,
    limit: usize
}

impl DensityQuery for AnywhereButInSetConstraintQuery {
    fn next_row(&mut self) -> Option<usize> {
        while self.current_row < self.limit {
            let t = self.current_row;

            self.current_row += 1;

            if self.not_in_rows.contains(&t) {
                continue
            } else {
                return Some(t);
            }
        }

        None
    }
}

impl AnywhereButInSetConstraint {
    fn eliminated_roots<F: PrimeField>(&self, column_domain: &Domain<F>) -> Vec<F> {
        let roots_generator = column_domain.generator;

        let mut rows = self.not_in_rows.clone();
        rows.sort();

        let mut roots = vec![];
        for num in rows.iter() {
            let x = roots_generator.pow([*num as u64]);
            roots.push(x);               
        }

        roots
    }
}

impl<F: PrimeField> ConstraintDensity<F> for AnywhereButInSetConstraint {
    fn calculate_density_query(&self, num_rows_in_trace: usize) -> Box<dyn DensityQuery> {
        assert!(self.not_in_rows.len() > 0);
        let max = *self.not_in_rows.iter().max().unwrap();
        assert!(max < num_rows_in_trace);

        let mut rows = self.not_in_rows.clone();
        rows.sort();

        let query = AnywhereButInSetConstraintQuery {
            not_in_rows: rows,
            current_row: 0,
            limit: num_rows_in_trace
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
        let divisor_degree = self.not_in_rows.len();

        let roots = self.eliminated_roots(&column_domain);

        let roots_iter = roots.iter();

        let column_domain_generator = column_domain.generator;
        let evaluation_domain_generator = evaluation_domain.generator;
        let multiplicative_generator = F::multiplicative_generator();

        let evaluation_domain_size = evaluation_domain.size as usize;
        // assert!(evaluation_domain_size == precomputations.coset.len());

        // strategy: 
        // - pick evaluation domain (will depend on constraint power) and column domain
        // - evaluate polynomial {X}^T - 1 where T is a size of the column domain over the coset of evaluation domain
        // - now we need to multiply those values by ({X} - root) for roots of this constraint
        // - first inverse {X}^T - 1 values
        // - for each of these values evaluate ({X} - root) and muptiply by this value

        // TODO: check if we need precomputation at all

        let mut vanishing = multiplicative_generator.pow([column_domain.size]);
        vanishing.sub_assign(&F::one());

        let vanishish_inverse = vanishing.inverse().unwrap();

        let inverse_divisors = vec![vanishish_inverse; column_domain.size as usize];

        let mut inverse_divisors = Polynomial::<F, Values>::from_values(inverse_divisors)?;

        // remove roots
        worker.scope(inverse_divisors.size(), |scope, chunk| {
            for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
                let roots_iter_outer = roots_iter.clone();
                scope.spawn(move |_| {
                    let mut x = column_domain_generator.pow([(i*chunk) as u64]);
                    x.mul_assign(&multiplicative_generator);
                    for v in inv_divis.iter_mut() {
                        let mut d = *v;
                        for root in roots_iter_outer.clone() {
                            // (X - root)
                            let mut tmp = x;
                            tmp.sub_assign(&root);
                            d.mul_assign(&tmp);
                        } 
                        // 1 / ( (X^T-1) / (X - 1)(X - omega)(...) ) =  (X - 1)(X - omega)(...) / (X^T-1)
                        *v = d;

                        x.mul_assign(&column_domain_generator);
                    }
                });
            }
        });

        inverse_divisors.batch_inversion(&worker)?;

        let divisor = inverse_divisors.icoset_fft(&worker);

        let mut divisor_lde = divisor.coset_lde(&worker, (evaluation_domain.size / column_domain.size) as usize)?;
        
        divisor_lde.batch_inversion(&worker)?;

        Ok((divisor_lde, divisor_degree))
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
        
        let mut vanishing_for_column_domain = at.pow([column_domain.size]);
        vanishing_for_column_domain.sub_assign(&F::one());

        let mut result = vanishing_for_column_domain.inverse().unwrap();

        for root in self.eliminated_roots(&column_domain).iter() {
            // (X - root)
            let mut tmp = at;
            tmp.sub_assign(&root);
            result.mul_assign(&tmp);
        } 

        Ok(result)
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