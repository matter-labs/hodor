use ff::PrimeField;
use crate::air::*;
use super::InstanceProperties;

pub trait IntoARP<F: PrimeField> {
    // return full trace, trace constraints and boundary constraints
    // registers should be remapped to have uniform structure, so ARP knows nothing about
    // what is a meaning of particular register
    fn into_arp(self) -> (Option<Vec<Vec<F>>>, InstanceProperties<F>);
}

// // ARP works with remapped registers and no longer cares about their meaning
// pub struct IntoARPInstance<F: PrimeField> {
//     pub witness: Option<Vec<Vec<F>>>,
//     pub num_rows: usize,
//     pub num_registers: usize,
//     pub constraints: Vec<Constraint<F>>,
//     pub boundary_constraints: Vec<BoundaryConstraint<F>>
// }