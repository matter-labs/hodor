use ff::PrimeField;
use crate::SynthesisError;

pub struct Domain<F: PrimeField> {
    pub size: u64,
    pub power_of_two: u64,
    pub generator: F,
}

impl<F: PrimeField> Domain<F> {
    pub fn new_for_size(size: u64) -> Result<Self, SynthesisError> {
        let power_of_two = size.next_power_of_two();
        let max_power_of_two = F::S as u64;
        if power_of_two > max_power_of_two {
            return Err(SynthesisError::Error);
        }

        let mut generator = F::root_of_unity();
        for _ in power_of_two..max_power_of_two {
            generator.square()
        }

        let size = 1u64 << power_of_two;

        Ok(Self {
            size: size,
            power_of_two: power_of_two,
            generator: generator
        })
    }
}