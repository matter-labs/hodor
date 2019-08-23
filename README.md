# Hodor

Open source zkSTARKs implementataion. Initial focus should be on highly efficient FRI reduction, with later developement of tools for AIR representation synthesis.

## Features
- [x] Full feature set (ARP + ALI) for formulation with single witness per register
- [x] DEEL-ALI as baseline for efficiency
- [x] Multicore (including fancy FFT strategies)
- [x] Radix-4 FFT
- [x] ZK by no quering from the original domains, but LDEs only. Additional masking by non-constrained elements of witness applies for free
- [ ] Mixed-radix FFT (2 and 4) for arbitrary domain lengths
- [ ] Constant registers optimization (WIP)
- [ ] One-shot and sparse constraints implementation (WIP)
- [x] Prover and verifier with precomputations at initialization time
- [ ] Carefull use of precomputations (WIP)
- [ ] Serialization formats
- [ ] Check if IOP Merkle trees should commit to natural domain indexing in addition to the evaluation result itself 
- [ ] Proof size optimization with coset combining
- [ ] Sparse FRI with less commitments to intermediate values
- [ ] Single FRI polynomial composition (WIP), pure proof size optimization that does not affect correctness and reduced proof speed
- [ ] Gadget library (WIP)
- [ ] Gadget composition synthesis-time optimizer and DSL for it
- [ ] Perfect ZK by masking polynomials

## License

Dual MIT/Apache-2.0

## Authors

- Alex Vlasov, [shamatar](https://github.com/shamatar)
- Konstantin Panarin, [Konstantce](https://github.com/Konstantce)