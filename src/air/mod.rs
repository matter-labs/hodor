mod test_trace_system;

use ff::{
    PrimeField,
};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum Register {
    ProgramCounter(usize),
    Register(usize),
    Constant(usize),
    Aux(usize)
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum ConstraintDensity {
    Dense,
    Sparse(usize)
}

pub struct UnivariateConstraintTerm<F: PrimeField>(F, (Register, usize), u64); // coeff, register, power
pub struct PolyvariateConstraintTerm<F: PrimeField>(F, Vec<UnivariateConstraintTerm<F>>); // coeff, terms
pub struct PolynomialConstraint<F: PrimeField>(F, Vec<PolynomialConstraintTerm<F>>); // linera combination of terms

pub enum PolynomialConstraintTerm<F: PrimeField> {
    Univariate(UnivariateConstraintTerm<F>),
    Polyvariate(PolyvariateConstraintTerm<F>)
}

pub trait TraceSystem<F: PrimeField> {
    // adds a new register to the system
    fn allocate_register(
        &mut self, 
        name: String
    ) -> Result<Register, String>;
    // adds a new constant register
    fn allocate_constant_register<CWF>(
        &mut self, 
        name: String,
        f: CWF
    ) -> Result<Register, String> where CWF: 'static + FnOnce(usize) -> Result<(F, bool), ()>;
    // tries to get aux register
    fn get_aux_register(
        &mut self, 
        register: usize
    ) -> Result<Register, String>;
    // adds constraint and 
    fn add_constraint<WF>(
        &mut self, 
        step: usize, 
        constraint: PolynomialConstraintTerm<F>, 
        density: ConstraintDensity, 
        value: WF,
    ) -> Result<(), String> where WF: 'static + FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, ()>;
    fn add_boundary_constraint(
        &mut self, 
        name: String,
        register: Register, 
        step: usize
    ) -> Result<(), String>;

    fn step(&mut self, num_steps: usize) -> Result<(), String>;
    // fn finalize(&mut self) -> Result<(), String>;
}


pub trait IntoAIR {
    fn trace<F: PrimeField, TC: TraceSystem<F>>(self, tracer: &mut TC) -> Result<(), ()>;
}

