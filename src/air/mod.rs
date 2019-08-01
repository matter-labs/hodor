mod test_trace_system;
mod constraint;

pub use constraint::*;
pub use test_trace_system::*;

use ff::{
    PrimeField,
};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum Register {
    ProgramCounter(usize),
    Register(usize),
    Constant(usize),
    Aux(usize)
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct DenseConstraint {
    pub start_at: usize,
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct RepeatedConstraint {
    pub start_at: usize,
    pub interval: usize
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct SparseConstraint {
    pub rows: Vec<usize>
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum ConstraintDensity {
    Dense(DenseConstraint),
    Repeated(RepeatedConstraint),
    Sparse(SparseConstraint),
}

impl std::default::Default for ConstraintDensity {
    fn default() -> Self {
        ConstraintDensity::Dense(DenseConstraint::default())
    }
}

impl DenseConstraint {
    pub fn new(starting_at: usize) -> Self {
        DenseConstraint{
            start_at: starting_at
        }
    }
}

impl std::default::Default for DenseConstraint {
    fn default() -> Self {
        DenseConstraint{
            start_at: 0
        }
    }
}

impl RepeatedConstraint {
    pub fn new(starting_at: usize, interval: usize) -> Self {
        assert!(interval != 1 && interval != 0, "Repeated constraint with interval 1 should be DenseConstraint");
        
        RepeatedConstraint{
            start_at: starting_at,
            interval: interval
        }
    }
}

impl std::default::Default for RepeatedConstraint {
    fn default() -> Self {
        RepeatedConstraint{
            start_at: 0,
            interval: 2
        }
    }
}

impl SparseConstraint{
    pub fn new(at_rows: Vec<usize>) -> Self {
        assert!(at_rows.len() > 0, "Interval constraint constraint must hold somewhere");
        
        SparseConstraint{
            rows: at_rows
        }
    }
}

impl std::default::Default for SparseConstraint {
    fn default() -> Self {
        SparseConstraint{
            rows: vec![]
        }
    }
}

#[derive(Debug)]
pub enum TracingError {
    Error,
}

use std::fmt;
use std::error::Error;

impl Error for TracingError {
    fn description(&self) -> &str {
        match *self {
            TracingError::Error => "General error for now",
        }
    }
}

impl fmt::Display for TracingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "{}", self.description())
    }
}

pub trait TraceSystem<F: PrimeField> {
    // adds a new register to the system
    fn allocate_register(
        &mut self, 
        name: String
    ) -> Result<Register, TracingError>;
    fn get_register(
        &self,
        step: usize,
        register: Register
    ) -> Result<F, TracingError>;
    // adds a new constant register
    fn allocate_constant_register<CWF>(
        &mut self, 
        name: String,
        f: CWF
    ) -> Result<Register, TracingError> where CWF: 'static + FnOnce(usize) -> Result<(F, bool), TracingError>;
    // tries to get aux register
    fn allocate_aux_register(
        &mut self, 
    ) -> Result<Register, TracingError>;
    // adds constraint and 
    fn add_constraint<WF>(
        &mut self, 
        constraint: Constraint<F>, 
        value: WF,
    ) -> Result<(), TracingError> where WF: 'static + FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, TracingError>;

    fn add_constraint_with_witness<WF>(
        &mut self, 
        constraint: Constraint<F>, 
        value: WF,
    ) -> Result<(), TracingError> where WF: 'static + Fn(&Self) -> Result<Vec<(F, Register, usize)>, TracingError>;
    
    fn add_boundary_constraint(
        &mut self, 
        name: String,
        register: Register, 
        at_step: usize,
        value: Option<F>,
    ) -> Result<(), TracingError>;

    fn step(&mut self, num_steps: usize) -> Result<(), TracingError>;
    fn get_step_number(&self) -> Result<usize, TracingError>;
    // fn finalize(&mut self) -> Result<(), String>;
}


pub trait IntoAIR<F: PrimeField> {
    fn trace<TC: TraceSystem<F>>(self, tracer: &mut TC) -> Result<(), TracingError>;
}

