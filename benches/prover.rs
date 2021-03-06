use criterion::{AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use ff::{Field, PrimeField};
use hodor::arp::{InstanceProperties, IntoARP};
use hodor::experiments::Fr;
use hodor::iop::optimized_iop::cubic_vdf::CubicVDF;
use hodor::iop::optimized_iop::fri::{
    FRIProof, FRIProofPrototype, FriProof, FriProofPrototype, NaiveFriIop,
};
use hodor::iop::optimized_iop::iop::TrivialBlake2sIOP;
use hodor::iop::optimized_iop::prover::Prover;
use hodor::transcript::Blake2sTranscript;
use hodor::{arp::PerRegisterARP, domains::Domain};

use hodor::experiments::cubic_vdf::CubicVDF as OldCubicVdf;
use hodor::fft::cooley_tukey_ntt::*;
use hodor::fri::{
    FRIProof as OldFRIProof, FRIProofPrototype as OldFRIProofPrototype, FriProof as OldFriProof,
    FriProofPrototype as OldFriProofPrototype, NaiveFriIop as OldNaiveFriIop,
};
use hodor::iop::blake2s_trivial_iop::TrivialBlake2sIOP as OldTrivialBlake2sIOP;
use hodor::prover::Prover as OldProver;

fn init_old_prover<F: PrimeField>(
    num_operations: usize,
    lde_factor: usize,
    output_at_degree_plus_one: usize,
) -> (
    Vec<Vec<F>>,
    OldProver<
        F,
        Blake2sTranscript<F>,
        OldTrivialBlake2sIOP<F>,
        OldFRIProofPrototype<F, OldTrivialBlake2sIOP<F>>,
        OldFRIProof<F, OldTrivialBlake2sIOP<F>>,
        OldNaiveFriIop<F, OldTrivialBlake2sIOP<F>>,
        PerRegisterARP,
    >,
) {
    let vdf_instance = CubicVDF::<F> {
        start_c0: F::one(),
        start_c1: F::one(),
        // num_operations: (1 << 20) - 1
        num_operations: num_operations,
    };
    let (witness_polys, props) = vdf_instance.into_arp();
    let witness_polys = witness_polys.expect("some witness");

    let prover_instance = OldProver::<
        F,
        Blake2sTranscript<F>,
        OldTrivialBlake2sIOP<F>,
        OldFRIProofPrototype<F, OldTrivialBlake2sIOP<F>>,
        OldFRIProof<F, OldTrivialBlake2sIOP<F>>,
        OldNaiveFriIop<F, OldTrivialBlake2sIOP<F>>,
        PerRegisterARP,
    >::new(props.clone(), lde_factor, output_at_degree_plus_one)
    .expect("must work");

    // (vdf_instance, prover_instance)
    (witness_polys, prover_instance)
}

fn init_optimized_prover<F: PrimeField>(
    num_operations: usize,
    lde_factor: usize,
    output_at_degree_plus_one: usize,
) -> (
    // CubicVDF<F>,
    Vec<Vec<F>>,
    Prover<
        F,
        Blake2sTranscript<F>,
        TrivialBlake2sIOP<F>,
        FRIProofPrototype<F, TrivialBlake2sIOP<F>>,
        FRIProof<F, TrivialBlake2sIOP<F>>,
        NaiveFriIop<F, TrivialBlake2sIOP<F>>,
        PerRegisterARP,
    >,
) {
    let vdf_instance = CubicVDF::<F> {
        start_c0: F::one(),
        start_c1: F::one(),
        // num_operations: (1 << 20) - 1
        num_operations: num_operations,
    };
    let (witness_polys, props) = vdf_instance.into_arp();
    let witness_polys = witness_polys.expect("some witness");

    let prover_instance = Prover::<
        F,
        Blake2sTranscript<F>,
        TrivialBlake2sIOP<F>,
        FRIProofPrototype<F, TrivialBlake2sIOP<F>>,
        FRIProof<F, TrivialBlake2sIOP<F>>,
        NaiveFriIop<F, TrivialBlake2sIOP<F>>,
        PerRegisterARP,
    >::new(props.clone(), lde_factor, output_at_degree_plus_one)
    .expect("must work");

    // (vdf_instance, prover_instance)
    (witness_polys, prover_instance)
}

fn init_prover_with_square_root_fft<F: PrimeField>(
    num_operations: usize,
    lde_factor: usize,
    output_at_degree_plus_one: usize,
) -> (
    Vec<Vec<F>>,
    Prover<
        F,
        Blake2sTranscript<F>,
        TrivialBlake2sIOP<F>,
        FRIProofPrototype<F, TrivialBlake2sIOP<F>>,
        FRIProof<F, TrivialBlake2sIOP<F>>,
        NaiveFriIop<F, TrivialBlake2sIOP<F>>,
        PerRegisterARP,
    >,
) {
    let log_n_half = num_operations.trailing_zeros()/2;
    println!("log n half {}", log_n_half);
    assert!(log_n_half &1 == 0, "num operations should be squared");

    let vdf_instance = CubicVDF::<F> {
        start_c0: F::one(),
        start_c1: F::one(),
        // num_operations: (1 << 20) - 1
        num_operations: num_operations,
    };
    let (witness_polys, props) = vdf_instance.into_arp();
    let witness_polys = witness_polys.expect("some witness");

    let prover_instance = Prover::<
        F,
        Blake2sTranscript<F>,
        TrivialBlake2sIOP<F>,
        FRIProofPrototype<F, TrivialBlake2sIOP<F>>,
        FRIProof<F, TrivialBlake2sIOP<F>>,
        NaiveFriIop<F, TrivialBlake2sIOP<F>>,
        PerRegisterARP,
    >::new(props.clone(), lde_factor, output_at_degree_plus_one)
    .expect("must work");

    // (vdf_instance, prover_instance)
    (witness_polys, prover_instance)
}

pub fn group(crit: &mut Criterion) {
    let lde_factor = 16;
    let output_at_degree_plus_one = 1;

    let mut group = crit.benchmark_group("Cubic VDF Prover");
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    let num_operations = (8..24).step_by(2).map(|e| (1 << e)-1);

    for num_operation in num_operations {
        let (new_witness_polys, new_prover_instance) =
            init_prover_with_square_root_fft::<Fr>(num_operation, lde_factor, output_at_degree_plus_one);

        let size = new_witness_polys[0].len();
        let new_size = lde_factor*size.next_power_of_two();

        

        let log_size = new_size.trailing_zeros();
        let domain = Domain::<Fr>::new_for_size(new_size as u64).expect("a domain");
        let omega = domain.generator;
        let inner_size = 1 << (log_size / 2);
        
        let inner_generator = omega.pow(&[inner_size as u64]);
        assert_eq!(inner_generator.pow(&[inner_size as u64]), Fr::one());

        let precomputed_twiddle_factors =
            hodor::fft::strided_fft::utils::precompute_twiddle_factors(
                &inner_generator,
                inner_size as usize,
            );

        println!("size {} lde factor {} new size {} inner size {}", size, lde_factor, new_size, inner_size);
        // println!("size {} lde factor {} inner size {}", size, lde_factor, inner_size);

        let bench_id = BenchmarkId::new("prover-with-square-root-lde", num_operation);
        group.bench_with_input(bench_id, &new_witness_polys, |b, input| {
            b.iter(|| {
                new_prover_instance
                    .prove_with_square_root_fft(input.to_vec(), &precomputed_twiddle_factors)
            })
        });

        // classic lde

        let (old_witness_polys, old_prover_instance) =
            init_old_prover::<Fr>(num_operation, lde_factor, output_at_degree_plus_one);

        let bench_id = BenchmarkId::new("prover-with-classic-lde", num_operation);
        group.bench_with_input(bench_id, &old_witness_polys, |b, input| {
            b.iter(|| old_prover_instance.prove(input.to_vec()))
        });

        // ntt lde

        let (optimized_witness_polys, optimized_prover_instance) =
            init_optimized_prover::<Fr>(num_operation, lde_factor, output_at_degree_plus_one);

        

        let bitreversed_omegas =
            BitReversedOmegas::<Fr>::new_for_domain_size(size.next_power_of_two());

        let bench_id = BenchmarkId::new("prover-with-ntt-lde", num_operation);
        group.bench_with_input(bench_id, &optimized_witness_polys, |b, input| {
            b.iter(|| optimized_prover_instance.prove_with_ntt(input.to_vec(), &bitreversed_omegas))
        });

        
    }
}
