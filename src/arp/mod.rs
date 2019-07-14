use crate::air::*;
use crate::air::TestTraceSystem;
use crate::polynomials::*;
use ff::PrimeField;

mod prime_field_arp;

pub use prime_field_arp::*;

// ARP works with remapped registers and no longer cares about their meaning
pub struct ARP<F:PrimeField> {
    pub witness: Option<Vec<Vec<F>>>,
    pub witness_poly: Option<WitnessPolynomial<F>>,
    pub num_steps: usize,
    pub num_registers: usize,
    pub constraints: Vec<Constraint<F>>,
    pub boundary_constraints: Vec<BoundaryConstraint<F>>
}

pub enum WitnessPolynomial<F: PrimeField> {
    Single(Polynomial<F, Coefficients>),
    PerRegister(Vec<Polynomial<F, Coefficients>>),
}

pub trait IntoARP<F: PrimeField> {
    // return full trace, trace constraints and boundary constraints
    // registers should be remapped to have uniform structure, so ARP knows nothing about
    // what is a meaning of particular register
    fn into_arp(self) -> ARP<F>;
}

impl<F: PrimeField> IntoARP<F> for TestTraceSystem<F> {
    fn into_arp(self) -> ARP<F> {
        let num_pc_registers = self.pc_registers.len();
        let num_registers = self.registers.len();
        let num_aux_registers = self.aux_registers.len();
        let num_constant_registers = self.constant_registers.len();

        let total_registers = num_pc_registers + num_registers + num_aux_registers + num_constant_registers;
        let num_steps = self.current_step;

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

        ARP::<F> {
            witness,
            witness_poly: None,
            num_steps,
            num_registers,
            constraints,
            boundary_constraints,
        }
    }
}