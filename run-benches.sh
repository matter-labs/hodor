#!/bin/sh
cargo build --release --benches
BENCH_FILE=$(find ./target/release/deps -name "benchmarks-*" | grep -v "*.d*" | tail -n1)
"./${BENCH_FILE}" --bench