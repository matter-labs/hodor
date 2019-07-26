use ff::PrimeField;

use crate::air::*;
use crate::polynomials::*;
use crate::arp::*;
use crate::*;
use crate::fft::multicore::Worker;

use super::ALI;

#[derive(Debug)]
pub struct DeepALI<F: PrimeField> {
    pub f_poly: Polynomial<F, Coefficients>,
    pub f_poly_lde_values: Option<Polynomial<F, Values>>,
    pub g_poly: Polynomial<F, Coefficients>,
    pub g_poly_lde_values: Option<Polynomial<F, Values>>,
    pub h_1_poly: Option<Polynomial<F, Values>>,
    pub h_2_poly: Option<Polynomial<F, Values>>,
    pub all_masks: Vec<StepDifference<F>>,
    pub z_m: Option<Vec<F>>,
    pub f_at_z_m: Option<Vec<F>>
}

impl<F: PrimeField> From<ALI<F>> for DeepALI<F> {
    fn from(mut ali: ALI<F>) -> DeepALI<F> {
        let mut all_masks = vec![];
        for k in ali.mask_applied_polynomials.keys() {
            all_masks.push(k.clone());
        }

        let f_poly = match ali.f_poly {
            WitnessPolynomial::Single(p) => {
                p
            },
            _ => {
                unimplemented!();
            }
        };

        let g_poly = ali.g_poly.take().expect("is some");

        DeepALI::<F> {
            f_poly: f_poly,
            f_poly_lde_values: None,
            g_poly: g_poly,
            g_poly_lde_values: None,
            h_1_poly: None,
            h_2_poly: None,
            all_masks: all_masks,
            z_m: None,
            f_at_z_m: None,
        }
    }
}

impl<F: PrimeField> DeepALI<F> {
    pub fn make_deep(
        &mut self,
        f_poly_lde_values: Polynomial<F, Values>,
        g_poly_lde_values: Polynomial<F, Values>,
        z: F,
    ) -> Result<(), SynthesisError>
    {
        let worker = Worker::new();

        // before this prover has already sent evaluations of f(x) at domain of 
        // size ~ rho * degree(f) (domain D) and evaluations of g(x) at domain of size
        // ~rho * degree(g) = rho * max(degree(constraints)) * degree(f) (domain D')
        // we expect z to be field element outsize of any domain and don't make any check
        // now we need to evaluate Z = \prod (X- z*M_i) at domain D, as well as 
        // U = interpolate(pairs((z*M_i, f(z*M_i))) at domain D
        // unfortunately M_i are not forming any "good" structure, so we evaluate them
        // either naively by value (for Z), or make a Lagrange interpolation and evaluate (for U).
        // Then we construct h1 = (f(x) - U(x)) / Z(x) and h2 = (g(x) - g(z)) / (x-z) by values in 
        // domains D and D' reprectively

        // at the end of the day we need to evaluate LDE of g, 

        let z_m: Vec<_> = self.all_masks.iter().map(|m| {
            let mut tmp = z;
            match m {
                StepDifference::Mask(mask) => {
                    tmp.mul_assign(&mask);
                },
                _ => {
                    unreachable!();
                }
            }

            tmp
        }).collect();

        let f_at_z_m: Vec<_> = z_m.iter().map(|m| {
            self.f_poly.evaluate_at(&worker, *m)
        }).collect();

        let mut points = vec![];

        for (z_m, f_at_z_m) in z_m.iter()
                                .zip(f_at_z_m.iter()) 
        {
            points.push((*z_m, *f_at_z_m));
        }

        // now we can interpolate F(z_m) as U in DEEP-ALI notations
        let u_coeffs = crate::utils::poly::interpolate(&points[..]).expect("must exist");
        let u_poly_coeffs = Polynomial::from_coeffs(u_coeffs)?;

        let z_coeffs = Polynomial::from_roots(z_m.clone(), &worker)?;

        let d_size = f_poly_lde_values.size();
        let d_prime_size = g_poly_lde_values.size();

        // evaluate h1

        let mut h1_values = f_poly_lde_values;

        {
            let u_poly_size = u_poly_coeffs.size();
            let u_poly_values = u_poly_coeffs.lde(&worker, d_size / u_poly_size)?;
            h1_values.sub_assign(&worker, &u_poly_values);
        }

        {
            let z_poly_size = z_coeffs.size();
            let mut z_poly_values = z_coeffs.lde(&worker, d_size / z_poly_size)?;
            z_poly_values.batch_inversion(&worker)?;
            h1_values.mul_assign(&worker, &z_poly_values);
        }

        // evaluate h2

        let mut denominator = Polynomial::<F, Values>::from_values(vec![F::one(); d_prime_size])?;
        let g = denominator.omega;
        denominator.distribute_powers(&worker, g);
        worker.scope(denominator.size(), |scope, chunk| {
            for v in denominator.as_mut().chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        v.sub_assign(&z);
                    }
                });
            }
        });

        denominator.batch_inversion(&worker)?;

        let mut h2_values = g_poly_lde_values;
        let g_at_z = self.g_poly.evaluate_at(&worker, z);

        worker.scope(h2_values.size(), |scope, chunk| {
            for v in h2_values.as_mut().chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v.iter_mut() {
                        v.sub_assign(&g_at_z);
                    }
                });
            }
        });

        h2_values.mul_assign(&worker, &denominator);

        self.z_m = Some(z_m);
        self.f_at_z_m = Some(f_at_z_m);
        self.h_1_poly = Some(h1_values);
        self.h_2_poly = Some(h2_values);

        Ok(())
    }
}


#[test]
fn test_fib_conversion_into_deep_ali() {
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::arp::IntoARP;
    use crate::ali::ALI;
    use crate::fft::multicore::Worker;
    use crate::ali::deep_ali::*;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
    let mut arp = ARP::<Fr>::new(test_tracer);
    arp.route_into_single_witness_poly().expect("must work");

    let mut ali = ALI::from(arp);
    let alpha = Fr::from_str("123").unwrap();
    ali.calculate_g(alpha).expect("must work");

    let mut deep_ali = DeepALI::from(ali);
    let z = Fr::from_str("62").unwrap();

    let lde_factor = 2;

    let worker = Worker::new();

    let f_lde_values = deep_ali.f_poly.clone().lde(&worker, lde_factor).unwrap();
    let f_coeffs = f_lde_values.clone().ifft(&worker);
    let f_values = deep_ali.f_poly.clone().fft(&worker);
    for i in 0..f_values.size() {
        assert!(f_values.as_ref()[i] == f_lde_values.as_ref()[i*lde_factor]);
    }
    for i in 0..deep_ali.f_poly.size() {
        assert!(deep_ali.f_poly.as_ref()[i] == f_coeffs.as_ref()[i]);
    }
    for i in deep_ali.f_poly.size()..f_coeffs.size() {
        assert!(f_coeffs.as_ref()[i].is_zero());
    }

    let g_lde_values = deep_ali.g_poly.clone().lde(&worker, lde_factor).unwrap();
    let g_coeffs = g_lde_values.clone().ifft(&worker);
    let g_values = deep_ali.g_poly.clone().fft(&worker);
    for i in 0..g_values.size() {
        assert!(g_values.as_ref()[i] == g_lde_values.as_ref()[i*lde_factor]);
    }
    for i in 0..deep_ali.g_poly.size() {
        assert!(deep_ali.g_poly.as_ref()[i] == g_coeffs.as_ref()[i]);
    }
    for i in deep_ali.g_poly.size()..g_coeffs.size() {
        assert!(g_coeffs.as_ref()[i].is_zero());
    }

    deep_ali.make_deep(f_lde_values, g_lde_values, z).expect("must work");

    // let h1_0 = deep_ali.h_1_poly.as_ref().unwrap().as_ref()[0].into_repr().as_ref()[0];
    // let h1_1 = deep_ali.h_1_poly.as_ref().unwrap().as_ref()[1].into_repr().as_ref()[0];

    // let h2_0 = deep_ali.h_2_poly.as_ref().unwrap().as_ref()[0].into_repr().as_ref()[0];
    // let h2_1 = deep_ali.h_2_poly.as_ref().unwrap().as_ref()[1].into_repr().as_ref()[0];

    // assert_eq!(h1_0, 119);
    // assert_eq!(h1_1, 87);
    // assert_eq!(h2_0, 99);
    // assert_eq!(h2_1, 226);

    let h1_values = deep_ali.h_1_poly.take().unwrap();
    let h1_coeffs = h1_values.ifft(&worker);
    println!("{:?}", h1_coeffs.as_ref());

    let h2_values = deep_ali.h_2_poly.take().unwrap();
    let h2_coeffs = h2_values.ifft(&worker);
    println!("{:?}", h2_coeffs.as_ref());

    // println!("H1 = {:?}", deep_ali.h_1_poly);
    // println!("H2 = {:?}", deep_ali.h_2_poly);
}