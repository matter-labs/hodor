use super::*;

/*

This is an attemt to create an abstraction for AIR. Trace system tries to call the witness generation
function at the very time step to generate the witness value. Unfortunately this does not allow to have the following benefits

- To have multicore witness generation for independent registers
- To have cases where it's much better to generate witness "backward in time"


For these reasons such abstraction was considered not too useful and gadget library is much better if
proper gadget composition model is found

*/

pub struct TestTraceSystem<F: PrimeField> {
    pub pc_registers: Vec<String>,
    pub registers: Vec<String>,
    pub constant_registers: Vec<String>,
    pub aux_registers: Vec<String>,

    pub pc_registers_witness: Vec<Vec<F>>,
    pub registers_witness: Vec<Vec<F>>,
    pub constant_registers_witness: Vec<Vec<F>>,
    pub aux_registers_witness: Vec<Vec<F>>,

    register_generators: Vec<Box<FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, TracingError> > >,
    witness_generators: Vec<Box<Fn(&Self) -> Result<Vec<(F, Register, usize)>, TracingError> > >,
    constant_register_generators: Vec<Box<FnOnce(usize) -> Result<(F, bool), TracingError> > >,

    pub constraints: Vec<Constraint<F>>,
    pub boundary_constraints: Vec<BoundaryConstraint<F>>,
    pub current_step: usize
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
    fn allocate_aux_register(
        &mut self
    ) -> Result<Register, TracingError>
    {
        let num_registers = self.aux_registers.len();
        self.aux_registers.push(format!("Aux({})", num_registers));
        self.aux_registers_witness.push(vec![]);

        Ok(Register::Aux(num_registers))
    }

    fn add_constraint<WF>(
        &mut self, 
        constraint: Constraint<F>, 
        value: WF,
    ) -> Result<(), TracingError> where WF: 'static + FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, TracingError>
    {
        self.constraints.push(constraint);
        self.register_generators.push(Box::new(value));

        Ok(())
    }

    fn add_constraint_with_witness<WF>(
        &mut self, 
        constraint: Constraint<F>, 
        value: WF,
    ) -> Result<(), TracingError> where WF: 'static + Fn(&Self) -> Result<Vec<(F, Register, usize)>, TracingError> {
        self.constraints.push(constraint);
        self.witness_generators.push(Box::new(value));

        Ok(())
    }
    
    fn add_boundary_constraint(
        &mut self, 
        _name: String,
        register: Register, 
        at_step: usize,
        value: Option<F>
    ) -> Result<(), TracingError>
    {
        let constraint = BoundaryConstraint::<F> {
            register: register,
            at_row: at_step,
            value: value
        };
        self.boundary_constraints.push(constraint);

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


pub(crate) struct Fibonacci<F: PrimeField> {
    pub(crate) final_b: Option<u64>,
    pub(crate) at_step: Option<usize>,
    pub(crate) _marker: std::marker::PhantomData<F>
}

impl<F: PrimeField> IntoAIR<F> for Fibonacci<F> {
    fn trace<TC: TraceSystem<F>>(self, tracer: &mut TC) -> Result<(), TracingError> {
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

        let mut fib_constraint_0 = Constraint::default();

        let mut fib_constraint_1 = Constraint::default();

        let a_register_now = UnivariateTerm::<F> {
            coeff: F::one(),
            register: a_register,
            steps_difference: StepDifference::Steps(0),
            power: 1
        };

        let b_register_now = UnivariateTerm::<F> {
            coeff: F::one(),
            register: b_register,
            steps_difference: StepDifference::Steps(0),
            power: 1
        };

        let a_next_step = UnivariateTerm::<F> {
            coeff: F::one(),
            register: a_register,
            steps_difference: StepDifference::Steps(1),
            power: 1
        };

        let b_next_step = UnivariateTerm::<F> {
            coeff: F::one(),
            register: b_register,
            steps_difference: StepDifference::Steps(1),
            power: 1
        };

        /*
            |   T   |   A   |   B   |   C   |
            |   0   |   0   |   1   |   1   |
            |   1   |   1   |   2   |   3   |
            |   2   |   2   |   3   |   5   |
            |   3   |   3   |   5   |   8   |
            
            We need to constraints for Transition Function

            1. constraint 0:
            b_{i} = a{i+1} 
            b_{i} - a{i+1}  = 0

            2. constraint 1:
            a_{i} + b_{i} = b_{i+1}
            a_{i} + b_{i} - b_{i+1} = 0
        */

        // constraint for registed a_new = b;
        fib_constraint_0 -= b_register_now.clone(); 
        fib_constraint_0 += a_next_step;
        // b_new = a + b;
        
        // 
        fib_constraint_1 -= a_register_now; 
        fib_constraint_1 -= b_register_now; 
        fib_constraint_1 += b_next_step;

        // witness deriveation function gives witness value from corresponding trace step
        // tracer implicityl tracks current step.

        tracer.add_constraint_with_witness(fib_constraint_0, witness_derivation_function_0)?;
        tracer.add_constraint_with_witness(fib_constraint_1, witness_derivation_function_1)?;


        if self.final_b.is_some() {
            /*
                We also need three boundary constraints for initial and last values

                1. boundary_constaint_0 => initial value for register A
                2. boundary_constaint_1 => initial value for register B
                3. boundary_constaint_2 => last value for register B
            */
            let final_b = self.final_b.unwrap();
            let at_step = self.at_step.unwrap();

            let initial_a = F::one();
            let initial_b = F::one();
            let final_b = F::from_str(&final_b.to_string()).unwrap();

            tracer.add_boundary_constraint("Initial A".to_string(), a_register, 0, Some(initial_a))?;
            tracer.add_boundary_constraint("Initial B".to_string(), b_register, 0, Some(initial_b))?;
            tracer.add_boundary_constraint("Final B".to_string(), b_register, at_step, Some(final_b))?;
        }

        Ok(())
    }
}

impl<F: PrimeField> TestTraceSystem<F> {
    pub(crate) fn new() -> Self {
        Self {
            pc_registers: vec![],
            registers: vec![],
            constant_registers: vec![],
            aux_registers: vec![],
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

    pub(crate) fn calculate_witness(&mut self, a: u64, b: u64, steps: usize) {
        let a0 = F::from_str("1").unwrap();
        let b0 = F::from_str("1").unwrap();

        self.registers_witness[0].push(a0);
        self.registers_witness[1].push(b0);

        for i in 0..steps {
            println!("Step {}", i);
            for generator in self.witness_generators.iter() {
                let values = generator(&self).unwrap();
                for (value, register, step_delta) in values.into_iter() {
                    println!("Value of {:?} register after step {} (at row {}) = {}", register, i, i + step_delta, value);
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
fn test_fibonacci_air() {
    use crate::Fr;

    let fib = Fibonacci::<Fr> {
        final_b: Some(5),
        at_step: Some(3),
        _marker: std::marker::PhantomData
    };

    let mut test_tracer = TestTraceSystem::<Fr>::new();
    fib.trace(&mut test_tracer).expect("should work");
    test_tracer.calculate_witness(1, 1, 3);
}