extern crate criterion;
use criterion::{Criterion, Bencher, PlotConfiguration, AxisScale, Throughput, BenchmarkId};
use std::convert::*;

pub fn bench_for_parameter_set<F, T: Copy + Clone + TryInto<u64> + 'static>(c: &mut Criterion, id: &str, with_throughput: bool, parameters: &[T], mut f: F)
where
    F: FnMut(&mut Bencher<'_>, &T),
    <T as TryInto<u64>>::Error: std::fmt::Debug
{
    let mut group = c.benchmark_group(id);
    for param in parameters.iter() {
        let param_as_u64: u64 = param.clone().try_into().unwrap_or_else(|e| panic!(format!("invalid parameter, can not convert to u64: {:?}", e)));
        if with_throughput {
            group.throughput(Throughput::Elements(param_as_u64));
        }
        // let pretty_description = format!("Parameter = {}", param_as_u64);
        group.bench_with_input(BenchmarkId::from_parameter(param_as_u64), param, |b, param| {
            f(b, param)
        });
    }
    group.finish();
}

pub fn bench_for_log_parameter_set<F, T: Copy + Clone + TryInto<u64> + 'static>(c: &mut Criterion, id: &str, with_throughput: bool, parameter_logs: &[T], mut f: F)
where
    F: FnMut(&mut Bencher<'_>, &T),
    <T as TryInto<u64>>::Error: std::fmt::Debug
{
    let mut group = c.benchmark_group(id);
    let plot_config = PlotConfiguration::default()
        .summary_scale(AxisScale::Logarithmic);
    group.plot_config(plot_config);
    for param in parameter_logs.iter() {
        let log_param_as_u64: u64 = param.clone().try_into().unwrap_or_else(|e| panic!(format!("invalid parameter, can not convert to u64: {:?}", e)));
        let param_as_u64: u64 = 1 << log_param_as_u64;
        if with_throughput {
            group.throughput(Throughput::Elements(param_as_u64));
        }
        // let pretty_description = format!("Parameter = {}, log2 of parametter = {}", param_as_u64, log_param_as_u64);
        let pretty_description = format!("2^{} = {}", log_param_as_u64, param_as_u64);
        group.bench_with_input(BenchmarkId::from_parameter(pretty_description), &param, |b, param| {
            f(b, param)
        });
    }
    group.finish();
}

pub fn log_parametrized_comparison_benchmark<T: Copy + Clone + TryInto<u64> + 'static, U, P>(
    c: &mut Criterion, 
    id: &str, 
    with_throughput: bool, 
    parameter_logs: &[T], 
    preparation_function: P,
    mut benchmarking_functions: Vec<(String, Box<dyn FnMut(&mut Bencher<'_>, &U) + 'static>)>,
) where
    P: Fn(&T) -> U,
    <T as TryInto<u64>>::Error: std::fmt::Debug
{
    let mut group = c.benchmark_group(id);
    let plot_config = PlotConfiguration::default()
        .summary_scale(AxisScale::Logarithmic);
    group.plot_config(plot_config);
    for param in parameter_logs.iter() {
        let log_param_as_u64: u64 = param.clone().try_into().unwrap_or_else(|e| panic!(format!("invalid parameter, can not convert to u64: {:?}", e)));
        let param_as_u64: u64 = 1 << log_param_as_u64;
        if with_throughput {
            group.throughput(Throughput::Elements(param_as_u64));
        }
        let parameters = preparation_function(param);
        for (name, executor) in benchmarking_functions.iter_mut() {
            let pretty_description = format!("2^{} = {}", log_param_as_u64, param_as_u64);
            let id = BenchmarkId::new(&*name, pretty_description);
            group.bench_with_input(id, &parameters, |b, param| {
                executor(b, param)
            });
        }

    }
    group.finish();
}