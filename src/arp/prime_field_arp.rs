use crate::air::{Register, PolynomialConstraint, ConstraintDensity, PolynomialConstraintTerm, PolyvariateConstraintTerm, UnivariateConstraintTerm};
use ff::PrimeField;
use super::IntoARP;

// ARP works with remapped registers and no longer cares about their meaning
pub struct ARP<F: PrimeField> {
    pub witness: Vec<Vec<F>>,
    pub constraints: Vec<(usize, PolynomialConstraint<F>, ConstraintDensity)>,
    pub boundaty_constraints: Vec<(Register, usize, F)>
}


impl<F: PrimeField> ARP<F> {
    pub(crate) fn new<I: IntoARP<F>>(f: I) -> Self {
        let (witness, constraints, boundaty_constraints) = f.into_arp();

        Self {
            witness,
            constraints,
            boundaty_constraints
        }
    }
}

#[test]
fn test_fib_conversion() {
    use crate::Fr;
    use crate::air::Fibonacci;
    use crate::air::TestTraceSystem;
    use crate::air::IntoAIR;
    use crate::arp::IntoARP;

    let fib = Fibonacci::<Fr> {
        first_a: Some(1),
        first_b: Some(1),
        final_a: Some(3),
        at_step: Some(2),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 5);
    let arp = ARP::<Fr>::new(test_tracer);
}

