use ff::{Field};


// multiply coefficients of the polynomial by the scalar
pub fn mul_polynomial_by_scalar<F: Field>(a: &mut [F], b: F) {
    use crate::fft::multicore::Worker;

    let worker = Worker::new();

    worker.scope(a.len(), |scope, chunk| {
        for a in a.chunks_mut(chunk)
        {
            scope.spawn(move |_| {
                for a in a.iter_mut() {
                    a.mul_assign(&b);
                }
            });
        }
    });
}

// elementwise add coeffs of one polynomial with coeffs of other, that are 
// first multiplied by a scalar 
pub fn mul_add_polynomials<F: Field>(a: &mut [F], b: &[F], c: F) {
    use crate::fft::multicore::Worker;

    let worker = Worker::new();

    assert_eq!(a.len(), b.len());

    worker.scope(a.len(), |scope, chunk| {
        for (a, b) in a.chunks_mut(chunk).zip(b.chunks(chunk))
        {
            scope.spawn(move |_| {
                for (a, b) in a.iter_mut().zip(b.iter()) {
                    let mut r = *b;
                    r.mul_assign(&c);

                    a.add_assign(&r);
                }
            });
        }
    });
}

extern crate crossbeam;
use self::crossbeam::channel::{unbounded};

pub fn evaluate_at_consequitive_powers<'a, F: Field> (
    coeffs: &[F],
    first_power: F,
    base: F
) -> F
    {
    use crate::fft::multicore::Worker;

    let (s, r) = unbounded();

    let worker = Worker::new();

    worker.scope(coeffs.len(), |scope, chunk| {
        for (i, coeffs) in coeffs.chunks(chunk).enumerate()
        {
            let s = s.clone();
            scope.spawn(move |_| {
                let mut current_power = base.pow(&[(i*chunk) as u64]);
                current_power.mul_assign(&first_power);

                let mut acc = F::zero();

                for p in coeffs {
                    let mut tmp = *p;
                    tmp.mul_assign(&current_power);
                    acc.add_assign(&tmp);

                    current_power.mul_assign(&base);
                }

                s.send(acc).expect("must send");
            });
        }
    });

    drop(s);

    // all threads in a scope have done working, so we can safely read
    let mut result = F::zero();

    loop {
        if r.is_empty() {
            break;
        }
        let value = r.recv().expect("must not be empty");
        result.add_assign(&value);
    }

    result
}

/// Perform a Lagrange interpolation for a set of points
/// It's O(n^2) operations, so use with caution
pub fn interpolate<F: Field>(
    points: &[(F, F)]
) -> Option<Vec<F>> {
    let max_degree_plus_one = points.len();
    assert!(max_degree_plus_one >= 2, "should interpolate for degree >= 1");
    let external_iter = points.clone().into_iter();
    let internal = points.clone();
    let mut coeffs = vec![F::zero(); max_degree_plus_one];
    for (k, p_k) in external_iter.enumerate() {
        let (x_k, y_k) = p_k;
        // coeffs from 0 to max_degree - 1
        let mut contribution = vec![F::zero(); max_degree_plus_one];
        let mut demoninator = F::one();
        let mut max_contribution_degree = 0;
        for (j, p_j) in internal.iter().enumerate() {
            let (x_j, _) = p_j;
            if j == k {
                continue;
            }

            let mut diff = x_k.clone();
            diff.sub_assign(&x_j);
            demoninator.mul_assign(&diff);

            if max_contribution_degree == 0 {
                max_contribution_degree = 1;
                contribution.get_mut(0).expect("must have enough coefficients").sub_assign(&x_j);
                contribution.get_mut(1).expect("must have enough coefficients").add_assign(&F::one());
            } else {
                let mul_by_minus_x_j: Vec<F> = contribution.iter().map(|el| {
                    let mut tmp = el.clone();
                    tmp.mul_assign(&x_j);
                    tmp.negate();

                    tmp
                }).collect();

                contribution.insert(0, F::zero());
                contribution.truncate(max_degree_plus_one);

                assert_eq!(mul_by_minus_x_j.len(), max_degree_plus_one);
                for (i, c) in contribution.iter_mut().enumerate() {
                    let other = mul_by_minus_x_j.get(i).expect("should have enough elements");
                    c.add_assign(&other);
                }
            }
        }

        demoninator = demoninator.inverse().expect("denominator must be non-zero");
        for (i, this_contribution) in contribution.into_iter().enumerate() {
            let c = coeffs.get_mut(i).expect("should have enough coefficients");
            let mut tmp = this_contribution;
            tmp.mul_assign(&demoninator);
            tmp.mul_assign(&y_k);
            c.add_assign(&tmp);
        }

    }

    Some(coeffs)
}