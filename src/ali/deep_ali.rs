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
    pub all_masks: Vec<StepDifference<F>>
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
            all_masks: all_masks
        }
    }
}

impl<F: PrimeField> DeepALI<F> {
    pub fn make_deep(
        &mut self,
        lde_factor: usize,
        z: F
    ) -> Result<(), SynthesisError>
    {
        let worker = Worker::new();

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

        for (z_m, f_at_z_m) in z_m.iter()
                                .zip(f_at_z_m.iter()) 
        {
            println!("F at zM: zM = {}, F = {}", z_m, f_at_z_m);
        }

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

    deep_ali.make_deep(16, z).expect("must work");
}