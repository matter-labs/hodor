use crate::air::{Register, PolynomialConstraint, ConstraintDensity, PolynomialConstraintTerm, PolyvariateConstraintTerm, UnivariateConstraintTerm};
use crate::air::TestTraceSystem;
use ff::PrimeField;

mod prime_field_arp;

pub use prime_field_arp::*;

pub trait IntoARP<F: PrimeField> {
    // return full trace, trace constraints and boundary constraints
    // registers should be remapped to have uniform structure, so ARP knows nothing about
    // what is a meaning of particular register
    fn into_arp(self) -> (Vec<Vec<F>>, Vec<(usize, PolynomialConstraint<F>, ConstraintDensity)>, Vec<(Register, usize, F)>);
}

impl<F: PrimeField> IntoARP<F> for TestTraceSystem<F> {
    fn into_arp(self) -> (Vec<Vec<F>>, Vec<(usize, PolynomialConstraint<F>, ConstraintDensity)>, Vec<(Register, usize, F)>) {
        let num_pc_registers = self.pc_registers.len();
        let num_registers = self.registers.len();
        let num_aux_registers = self.aux_registers.len();
        let num_constant_registers = self.constant_registers.len();

        let register_remap_constant = num_pc_registers;
        let aux_register_remap_constant = register_remap_constant + num_registers;
        let constant_register_remap_constant = aux_register_remap_constant + num_aux_registers;

        let mut witness = vec![];

        fn remap(
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
            term: UnivariateConstraintTerm<F>,
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) -> UnivariateConstraintTerm<F> {
            let UnivariateConstraintTerm(c, reg, power) = term;
            let (reg, time_step) = reg;
            let reg = remap(reg, register_remap_constant, aux_register_remap_constant, constant_register_remap_constant);
            UnivariateConstraintTerm(c, (reg, time_step), power)
        }

        fn remap_term<F: PrimeField>(
            term: PolynomialConstraintTerm<F>,
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) -> PolynomialConstraintTerm<F> {
            let remapped = match term {
                PolynomialConstraintTerm::Univariate(uni_term) => {
                    PolynomialConstraintTerm::Univariate(
                        remap_univariate_term(
                            uni_term, 
                            register_remap_constant, 
                            aux_register_remap_constant, 
                            constant_register_remap_constant
                        )
                    )
                },
                PolynomialConstraintTerm::Polyvariate(poly_term) => {
                    let PolyvariateConstraintTerm(coeff, terms) = poly_term;
                    let terms: Vec<_> = terms.into_iter().map(|t| 
                        remap_univariate_term(t, register_remap_constant, aux_register_remap_constant, constant_register_remap_constant)
                    ).collect();

                    PolynomialConstraintTerm::Polyvariate(
                        PolyvariateConstraintTerm(coeff, terms)
                    )
                }
            };

            remapped
        }

        fn remap_constraint<F: PrimeField>(
            constraint: PolynomialConstraint<F>,
            register_remap_constant: usize, 
            aux_register_remap_constant: usize, 
            constant_register_remap_constant: usize
        ) -> PolynomialConstraint<F> {
            let PolynomialConstraint(coeff, terms) = constraint;
            let terms: Vec<_> = terms.into_iter().map(|t| 
                remap_term(t, register_remap_constant, aux_register_remap_constant, constant_register_remap_constant)
            ).collect();

            PolynomialConstraint(coeff, terms)
        }

        let constraints: Vec<_> = self.constraints.into_iter().map(|(step, c, density)| 
        {
            let c = remap_constraint(c, register_remap_constant, aux_register_remap_constant, constant_register_remap_constant);

            (step, c, density)
        }).collect();

        let boundary_constraints: Vec<_> = self.boundary_constraints.into_iter().map(|(reg, step, value)| 
        {
            let reg = remap(reg, register_remap_constant, aux_register_remap_constant, constant_register_remap_constant);

            (reg, step, value)
        }).collect();

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

        assert_eq!(witness.len(), num_pc_registers + num_registers + num_aux_registers + num_constant_registers);

        (witness, constraints, boundary_constraints)
    }
}