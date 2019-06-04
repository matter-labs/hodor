use super::*;

use std::ops::{AddAssign, SubAssign, MulAssign};

pub struct TestTraceSystem<F: PrimeField> {
    pc_registers: Vec<String>,
    registers: Vec<String>,
    constant_registers: Vec<String>,
    aux_registers: usize,
    pc_registers_witness: Vec<Vec<F>>,
    registers_witness: Vec<Vec<F>>,
    constant_registers_witness: Vec<Vec<F>>,
    register_generators: Vec<Box<FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, TracingError> > >,
    witness_generators: Vec<Box<Fn(&Self) -> Result<Vec<(F, Register, usize)>, TracingError> > >,
    constant_register_generators: Vec<Box<FnOnce(usize) -> Result<(F, bool), TracingError> > >,
    aux_registers_witness: Vec<Vec<F>>,
    constraints: Vec<(usize, PolynomialConstraint<F>, ConstraintDensity)>,
    boundary_constraints: Vec<(usize, Register)>,
    current_step: usize
}



impl<F: PrimeField> TraceSystem<F> for TestTraceSystem<F> {
// adds a new register to the system
    fn allocate_register(
        &mut self, 
        name: String,
    ) -> Result<Register, TracingError> 
    {
        let num_registers = self.registers.len();
        self.registers.push(name);
        self.registers_witness.push(vec![]);

        Ok(Register::Register(num_registers))
    }

    fn get_register(
        &self,
        step: usize,
        register: Register
    ) -> Result<F, TracingError> {
        match register {
            Register::Register(register_num) => {
                if register_num >= self.registers_witness.len() {
                    return Err(TracingError::Error);
                }
                let witness = &self.registers_witness[register_num];
                if step >= witness.len() {
                    return Err(TracingError::Error);
                }

                return Ok(witness[step]);
            },
            _ => {
                return Err(TracingError::Error);
            }
        }
    }
    // adds a new constant register
    fn allocate_constant_register<CWF>(
        &mut self, 
        name: String,
        f: CWF
    ) -> Result<Register, TracingError> where CWF: 'static + FnOnce(usize) -> Result<(F, bool), TracingError>
    {
        let num_registers = self.constant_registers.len();
        self.constant_registers.push(name);
        self.constant_registers_witness.push(vec![]);
        self.constant_register_generators.push(Box::new(f));

        Ok(Register::Constant(num_registers))
    }
    // tries to get aux register
    fn get_aux_register(
        &mut self, 
        register: usize
    ) -> Result<Register, TracingError>
    {
        if register >= self.aux_registers {
            self.aux_registers = register - 1;
            self.aux_registers_witness.resize(register, vec![]);
        }

        Ok(Register::Aux(register))
    }
    fn add_constraint<WF>(
        &mut self, 
        step: usize, 
        constraint: PolynomialConstraint<F>, 
        density: ConstraintDensity, 
        value: WF,
    ) -> Result<(), TracingError> where WF: 'static + FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, TracingError>
    {
        self.constraints.push((step, constraint, density));
        self.register_generators.push(Box::new(value));

        Ok(())
    }

    fn add_constraint_with_witness<WF>(
        &mut self, 
        step: usize, 
        constraint: PolynomialConstraint<F>, 
        density: ConstraintDensity, 
        value: WF,
    ) -> Result<(), TracingError> where WF: 'static + Fn(&Self) -> Result<Vec<(F, Register, usize)>, TracingError> {
        self.constraints.push((step, constraint, density));
        self.witness_generators.push(Box::new(value));

        Ok(())
    }
    
    fn add_boundary_constraint(
        &mut self, 
        _name: String,
        register: Register, 
        step: usize
    ) -> Result<(), TracingError>
    {
        self.boundary_constraints.push((step, register));

        Ok(())
    }

    fn step(&mut self, num_steps: usize) -> Result<(), TracingError> 
    {
        if num_steps == 0 {
            return Err(TracingError::Error);
        }
        self.current_step += num_steps;
        Ok(())
    }
    fn get_step_number(&self) -> Result<usize, TracingError> {
        Ok(self.current_step)
    }
    // fn finalize(&mut self) -> Result<(), String>;
} 


struct Fibonacci {
    first_a: Option<u64>,
    first_b: Option<u64>,
    final_a: Option<u64>,
    at_step: Option<usize>
}

impl IntoAIR for Fibonacci {
    fn trace<F: PrimeField, TC: TraceSystem<F>>(self, tracer: &mut TC) -> Result<(), TracingError> {
        let a_register = tracer.allocate_register("A".to_string())?;
        let b_register = tracer.allocate_register("B".to_string())?;

        let witness_derivation_function_0 = move |witness: &TC| -> Result<Vec<(F, Register, usize)>, TracingError> {
            let step_number = witness.get_step_number()?;
            let value = witness.get_register(step_number, b_register)?;
            let new_value = value;

            Ok(vec![(new_value, a_register, 1)])
        };

        let witness_derivation_function_1 = move |witness: &TC| -> Result<Vec<(F, Register, usize)>, TracingError> {
            let step_number = witness.get_step_number()?;
            let a_value = witness.get_register(step_number, a_register)?;
            let b_value = witness.get_register(step_number, b_register)?;
            let mut new_value = a_value;
            new_value.add_assign(&b_value);

            Ok(vec![(new_value, b_register, 1)])
        };

        let mut fib_constraint_0 = PolynomialConstraint::default();
        let mut fib_constraint_1 = PolynomialConstraint::default();
        let a_register_now = UnivariateConstraintTerm(F::one(), (a_register, 0), 1);
        let b_register_now = UnivariateConstraintTerm(F::one(), (b_register, 0), 1);
        let a_next_step = UnivariateConstraintTerm(F::one(), (a_register, 1), 1);
        let b_next_step = UnivariateConstraintTerm(F::one(), (b_register, 1), 1);
        // constraint for registed a_new = b;
        fib_constraint_0 += b_register_now.clone(); 
        fib_constraint_0 -= a_next_step;
        // b_new = a + b;
        fib_constraint_1 += a_register_now; 
        fib_constraint_1 += b_register_now; 
        fib_constraint_1 -= b_next_step;

        tracer.add_constraint_with_witness(0, fib_constraint_0, ConstraintDensity::Dense, witness_derivation_function_0)?;
        tracer.add_constraint_with_witness(0, fib_constraint_1, ConstraintDensity::Dense, witness_derivation_function_1)?;


        if self.final_a.is_some() {
            let a = self.first_a.unwrap();
            let b = self.first_b.unwrap();
            let final_a = self.final_a.unwrap();
            let at_step = self.at_step.unwrap();

            tracer.add_boundary_constraint("Initial A".to_string(), a_register, 0);
            tracer.add_boundary_constraint("Initial B".to_string(), b_register, 0);
            tracer.add_boundary_constraint("Final A".to_string(), a_register, at_step);
        }

        Ok(())
    }
}

impl<F: PrimeField> TestTraceSystem<F> {
    fn new() -> Self {
        Self {
            pc_registers: vec![],
            registers: vec![],
            constant_registers: vec![],
            aux_registers: 0,
            pc_registers_witness: vec![vec![]],
            registers_witness: vec![vec![]],
            constant_registers_witness: vec![vec![]],
            register_generators: vec![],
            witness_generators: vec![],
            constant_register_generators: vec![],
            aux_registers_witness: vec![vec![]],
            constraints: vec![],
            boundary_constraints: vec![],
            current_step: 0
        }
    }
    fn calculate_witness(&mut self, a: u64, b: u64, steps: usize) {
        let a0 = F::from_str("1").unwrap();
        let b0 = F::from_str("1").unwrap();

        self.registers_witness[0].push(a0);
        self.registers_witness[1].push(b0);

        for i in 0..steps {
            println!("Step {}", i);
            for generator in self.witness_generators.iter() {
                let values = generator(&self).unwrap();
                for (value, register, step_delta) in values.into_iter() {
                    println!("Value of {:?} register at step {} is {}", register, i + step_delta, value);
                    match register {
                        Register::Register(reg_num) => {
                            let reg_witness = &mut self.registers_witness[reg_num];
                            let at_step = i + step_delta;
                            if reg_witness.len() <= at_step {
                                reg_witness.resize(at_step + 1, F::zero());
                            }
                            reg_witness[at_step] = value;
                        },
                        _ => {
                            unreachable!();
                        }
                    }
                }
            }
            self.current_step += 1;
        }
    }
}

#[test]
fn test_fib() {
    use crate::Fr;

    let mut fib = Fibonacci {
        first_a: Some(1),
        first_b: Some(1),
        final_a: Some(3),
        at_step: Some(2)
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer);
    test_tracer.calculate_witness(1, 1, 5);
}