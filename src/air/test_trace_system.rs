use super::*;

pub struct TestTraceSystem<F: PrimeField> {
    pc_registers: Vec<String>,
    registers: Vec<String>,
    constant_registers: Vec<String>,
    aux_registers: usize,
    pc_registers_witness: Vec<Vec<F>>,
    registers_witness: Vec<Vec<F>>,
    constant_registers_witness: Vec<Vec<F>>,
    register_generators: Vec<Box<FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, ()> > >,
    constant_register_generators: Vec<Box<FnOnce(usize) -> Result<(F, bool), ()> > >,
    aux_registers_witness: Vec<Vec<F>>,
    constraints: Vec<(usize, PolynomialConstraint<F>, ConstraintDensity)>,
    // boundary_constraints: 
}

impl<F: PrimeField> TraceSystem<F> for TestTraceSystem<F> {
// adds a new register to the system
    fn allocate_register(
        &mut self, 
        name: String,
    ) -> Result<Register, String> 
    {
        let num_registers = self.registers.len();
        self.registers.push(name);
        self.registers_witness.push(vec![]);

        Ok(Register::Register(num_registers))
    }
    // adds a new constant register
    fn allocate_constant_register<CWF>(
        &mut self, 
        name: String,
        f: CWF
    ) -> Result<Register, String> where CWF: 'static + FnOnce(usize) -> Result<(F, bool), ()>
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
    ) -> Result<Register, String>
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
        constraint: PolynomialConstraintTerm<F>, 
        density: ConstraintDensity, 
        value: WF,
    ) -> Result<(), String> where WF: 'static + FnOnce(Vec<(F, Register, usize)>) -> Result<Vec<(F, Register, usize)>, ()>
    {
        Ok(())
    }


    fn add_boundary_constraint(
        &mut self, 
        name: String,
        register: Register, 
        step: usize
    ) -> Result<(), String>
    {
        Ok(())
    }

    fn step(&mut self, num_steps: usize) -> Result<(), String> 
    {
        Ok(())
    }
    // fn finalize(&mut self) -> Result<(), String>;
} 
