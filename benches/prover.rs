use criterion::{AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use ff::{Field, PrimeField};
use hodor::arp::PerRegisterARP;
use hodor::arp::{InstanceProperties, IntoARP};
use hodor::experiments::Fr;
use hodor::iop::optimized_iop::cubic_vdf::CubicVDF;
use hodor::iop::optimized_iop::fri::{
    FRIProof, FRIProofPrototype, FriProof, FriProofPrototype, NaiveFriIop,
};
use hodor::iop::optimized_iop::iop::TrivialBlake2sIOP;
use hodor::iop::optimized_iop::prover::Prover;
use hodor::transcript::Blake2sTranscript;

use hodor::experiments::cubic_vdf::CubicVDF as OldCubicVdf;
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

pub fn group(crit: &mut Criterion) {
    let lde_factor = 16;
    let output_at_degree_plus_one = 1;

    let mut group = crit.benchmark_group("Cubic VDF Prover");
    let new_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    group.plot_config(new_config);

    let num_operations = (8..16).map(|e| (1 << e) - 1);

    for num_operation in num_operations {
        let (old_witness_polys, old_prover_instance) =
            init_old_prover::<Fr>(num_operation, lde_factor, output_at_degree_plus_one);

        let bench_id = BenchmarkId::new("old-prover", num_operation);
        group.bench_with_input(bench_id, &old_witness_polys, |b, input| {
            b.iter(|| old_prover_instance.prove(input.to_vec()))
        });

        let (optimized_witness_polys, optimized_prover_instance) =
            init_optimized_prover::<Fr>(num_operation, lde_factor, output_at_degree_plus_one);

        let bench_id = BenchmarkId::new("optimized-prover", num_operation);
        group.bench_with_input(bench_id, &optimized_witness_polys, |b, input| {
            b.iter(|| optimized_prover_instance.prove(input.to_vec()))
        });
    }
}
