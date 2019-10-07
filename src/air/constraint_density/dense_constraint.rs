use super::*;

// Happens at every row including start_at
// Span is an auxilary information that will not allow to "wrap around" the trace
// and also can be used to stop earlier
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct DenseConstraint {
    pub start_at: usize,
    pub span: usize,
}


impl DenseConstraint {
    pub fn new(starting_at: usize, span: usize) -> Self {
        assert!(span >= 1, "Span >= 1");
        
        DenseConstraint{
            start_at: starting_at,
            span: span,
        }
    }

    pub fn new_wrapping(starting_at: usize, span: usize) -> Self {
        DenseConstraint{
            start_at: starting_at,
            span: span,
        }
    }
}

impl std::default::Default for DenseConstraint {
    fn default() -> Self {
        DenseConstraint{
            start_at: 0,
            span: 1,
        }
    }
}


pub(crate) struct DenseConstraintQuery {
    current_row: usize,
    limit: usize
}

impl DensityQuery for DenseConstraintQuery {
    fn next_row(&mut self) -> Option<usize> {
        if self.current_row < self.limit {
            let t = self.current_row;
            self.current_row += 1;

            Some(t)
        } else {
            None
        }
    }
}

impl DenseConstraint {
    fn eliminated_roots<F: PrimeField>(&self, column_domain: &Domain<F>) -> Vec<F> {
        let roots_generator = column_domain.generator;

        let start_at = self.start_at as u64;
        let span = self.span as u64;

        let mut roots = vec![];
        let mut root = F::one();
        for _ in 0..start_at {
            roots.push(root);
            root.mul_assign(&roots_generator);                
        }
        
        let last_step = column_domain.size - span;
        let mut root = roots_generator.pow([last_step]);
        for _ in last_step..column_domain.size {
            roots.push(root);
            root.mul_assign(&roots_generator);
        }

        roots
    }
}

impl<F: PrimeField> ConstraintDensity<F> for DenseConstraint {
    fn calculate_density_query(&self, num_rows_in_trace: usize) -> Box<dyn DensityQuery> {
        let start_at = self.start_at;
        let span = self.span;

        let limit = num_rows_in_trace - span;

        let query = DenseConstraintQuery {
            current_row: start_at,
            limit: limit
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
        let start_at = self.start_at;
        let span = self.span as u64;
        let mut divisor_degree = column_domain.size as usize;
        let divisor_domain_size = column_domain.size as usize;
        divisor_degree -= start_at as usize;
        divisor_degree -= (divisor_domain_size as usize) - num_rows_in_trace;
        divisor_degree -= span as usize;

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

        // // these are values at the coset
        // let mut inverse_divisors = Polynomial::<F, Values>::new_for_size(evaluation_domain_size)?;
        
        // // prepare for batch inversion
        // worker.scope(inverse_divisors.size(), |scope, chunk| {
        //     for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
        //         scope.spawn(move |_| {
        //             let mut x = evaluation_domain_generator.pow([(i*chunk) as u64]);
        //             x.mul_assign(&multiplicative_generator);
        //             for v in inv_divis.iter_mut() {
        //                 *v = x.pow([divisor_domain_size as u64]);
        //                 v.sub_assign(&F::one());

        //                 x.mul_assign(&evaluation_domain_generator);
        //             }
        //         });
        //     }
        // });

        // // now polynomial is filled with X^T - 1, and need to be inversed

        // inverse_divisors.batch_inversion(&worker)?;

        // // now do the evaluation

        // worker.scope(inverse_divisors.size(), |scope, chunk| {
        //     for (i, inv_divis) in inverse_divisors.as_mut().chunks_mut(chunk).enumerate() {
        //         let roots_iter_outer = roots_iter.clone();
        //         scope.spawn(move |_| {
        //             let mut x = evaluation_domain_generator.pow([(i*chunk) as u64]);
        //             x.mul_assign(&multiplicative_generator);
        //             for v in inv_divis.iter_mut() {
        //                 let mut d = *v;
        //                 for root in roots_iter_outer.clone() {
        //                     // (X - root)
        //                     let mut tmp = x;
        //                     tmp.sub_assign(&root);
        //                     d.mul_assign(&tmp);
        //                 } 
        //                 // 1 / ( (X^T-1) / (X - 1)(X - omega)(...) ) =  (X - 1)(X - omega)(...) / (X^T-1)
        //                 *v = d;

        //                 x.mul_assign(&evaluation_domain_generator);
        //             }
        //         });
        //     }
        // });

        // Ok((inverse_divisors, divisor_degree))
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