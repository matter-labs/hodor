use crate::air::*;
use crate::air::TestTraceSystem;
use crate::polynomials::*;
use crate::fft::multicore::Worker;
use crate::transcript::Transcript;
use crate::*;

use ff::PrimeField;
mod into_arp;
mod single_witness;
mod per_register;

pub use into_arp::*;

pub trait ARPType: 
    Sized 
    + Copy 
    + Clone
    {}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum SingleWitnessARP { }

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum PerRegisterARP { }

impl ARPType for SingleWitnessARP {}
impl ARPType for PerRegisterARP {}

pub trait ARP<F :PrimeField>: 
    Sized 
    + Clone 
    + std::fmt::Debug 
    {
        fn from_instance(
            instance: IntoARPInstance<F>, 
            worker: &Worker
        ) -> Result<Self, SynthesisError>;
    }


#[derive(Debug, Clone)]
pub enum WitnessPolynomial<F: PrimeField> {
    Single(Polynomial<F, Coefficients>),
    PerRegister(Vec<Polynomial<F, Coefficients>>),
}

impl<F: PrimeField> IntoARP<F> for TestTraceSystem<F> {
    fn into_arp(self) -> IntoARPInstance<F> {
        let num_pc_registers = self.pc_registers.len();
        let num_registers = self.registers.len();
        let num_aux_registers = self.aux_registers.len();
        let num_constant_registers = self.constant_registers.len();

        let total_registers = num_pc_registers + num_registers + num_aux_registers + num_constant_registers;
        let num_rows = self.current_step+1;

        let register_remap_constant = num_pc_registers;
        let aux_register_remap_constant = register_remap_constant + num_registers;
        let constant_register_remap_constant = aux_register_remap_constant + num_aux_registers;

        let mut witness = vec![];

        fn remap_register(
            register: Register, 
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) -> Register {
            let remapped = match register {
                Register::ProgramCounter(i) => {
                    Register::Register(i)
                },
                Register::Register(i) => {
                    Register::Register(i+register_remap_constant)
                },
                Register::Aux(i) => {
                    Register::Register(i+aux_register_remap_constant)
                },
                Register::Constant(i) => {
                    Register::Register(i+constant_register_remap_constant)
                }
            };

            remapped
        }

        fn remap_univariate_term<F: PrimeField>(
            term: &mut UnivariateTerm<F>,
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) {
            term.register = remap_register(
                term.register,
                register_remap_constant,
                aux_register_remap_constant, 
                constant_register_remap_constant
            );
        }

        fn remap_term<F: PrimeField>(
            term: &mut ConstraintTerm<F>,
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) {
            match term {
                ConstraintTerm::Univariate(ref mut t) => {
                    remap_univariate_term(
                        t, 
                        register_remap_constant, 
                        aux_register_remap_constant, 
                        constant_register_remap_constant
                    );
                },
                ConstraintTerm::Polyvariate(ref mut poly_term) => {
                    for mut t in poly_term.terms.iter_mut() {
                        remap_univariate_term(
                            &mut t, 
                            register_remap_constant, 
                            aux_register_remap_constant, 
                            constant_register_remap_constant
                        );
                    }
                }
            }
        }

        fn remap_constraint<F: PrimeField>(
            constraint: &mut Constraint<F>,
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) {
            for mut t in constraint.terms.iter_mut() {
                remap_term(
                    &mut t, 
                    register_remap_constant, 
                    aux_register_remap_constant, 
                    constant_register_remap_constant);
            }
        }

        let mut constraints = self.constraints;

        for mut c in constraints.iter_mut() {
            remap_constraint(
                &mut c, 
                register_remap_constant, 
                aux_register_remap_constant, 
                constant_register_remap_constant
            );
        }   

        let mut boundary_constraints = self.boundary_constraints;
        for c in boundary_constraints.iter_mut() {
            c.register = remap_register(
                c.register, 
                register_remap_constant, 
                aux_register_remap_constant, 
                constant_register_remap_constant
            );
        }   

        for r in self.pc_registers_witness.into_iter() {
            if r.len() != 0 {
                witness.push(r);
            }
        }

        for r in self.registers_witness.into_iter() {
            if r.len() != 0 {
                witness.push(r);
            }
        }

        for r in self.aux_registers_witness.into_iter() {
            if r.len() != 0 {
                witness.push(r);
            }
        }

        for r in self.constant_registers_witness.into_iter() {
            if r.len() != 0 {
                witness.push(r);
            }
        }

        assert_eq!(witness.len(), total_registers);

        let witness = if witness.len() == 0 {
            None
        } else {
            Some(witness)
        };

        IntoARPInstance::<F> {
            witness,
            num_rows,
            num_registers,
            constraints,
            boundary_constraints,
        }
    }
}