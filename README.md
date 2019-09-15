# Hodor

Open source zkSTARKs implementataion over prime fields. Initial focus should be on highly efficient prover implementation, with later developement of tools for AIR representation synthesis.

## Release details

Matter Labs has applied to present an open-source zkSTARK prover on Devcon5 with intended opening date of early/mid September. Due to publication of our transparent and quantum secure [commitment scheme](https://eprint.iacr.org/2019/1020) we think that even an alpha version of our prover/framework is now of a separate interest and will allow end-users to build an intuition on zkSTARK, FRI and IOP.

## Challenges

Stark arithmetization - AIR - is much more difficult to design for than traditional R1CS or gate-based arithmetizations from Groth16, Sonic, Bulletproofs or PLONK due to lack of "memory" in a form of ability to always address a previously declared variable and being intrinsically non-suitable for one-off computations. As a result creation of a programming and (much more difficult) some form of gadget-composability approach is a non-trivial task and here the spirit of open source programming and community should help to find an optimal track. More material about the essence of Stark constraints and their "density" will be published soon.

## Current state of affairs

Current code is not yet cleaned and should not be used for anything but proofs of concept, but is close to production grade in a sense that it's designed from scratch by following the publications (it's not a port of some other code), mainly DEEP-FRI paper, and it multicore optimized from the day one. It's also commented so understand the workflow of the prover starting from the ARP step.

In the prover it's assumed that single round of FRI is enough to reach the target soundness error (e.g. 100 bits of security). It's not always the case, so implementation of simpler interface for FRI parametrization is a top priority (see below).

Also, only "dense" constraints are implemented for now - constraints that affect every next row of the AIR trace with may be few rows skipped at the benining and at the end. Other types of well-computable densities, e.g. one that happens at every 2nd, 4th, 8th, etc. row of the trace (due to implementation being over multiplicative subgroup of the size 2^k in a prime field for FFT purposes) and some other, are not implemented yet.

## Prioritized TODO (contributions welcome)

- [ ] Implement (and most likely rework the API) for the existing constraint densities
- [ ] Provide more parametrizable FRI (more rounds, skipping intermediate commitment steps)
- [ ] Provide some trivial AIR toolkit (e.g. step by step tracer and witness generator)
- [ ] Decide on the programming model for gadget compositions
- [ ] Abstract away more constraint "densities" (divisors for ALI step)
- [ ] Cleanup traits and public interfaces
- [ ] Add native Rust serialization (`serde` based) for constraint systems or at least constraint densities for non-trivial cases 
- [ ] Automate a workflow for "constant" registers (lookup tables for e.g. Pedersed hash) 
- [ ] Would be good, but not strictly: Ethereum verifier example or (hard way) synthesiser

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
- [ ] Perfect ZK by masking polynomials

## License

Dual MIT/Apache-2.0

## Authors

- Alex Vlasov, [shamatar](https://github.com/shamatar)
- Konstantin Panarin, [Konstantce](https://github.com/Konstantce)