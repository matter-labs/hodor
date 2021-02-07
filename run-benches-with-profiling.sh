#!/bin/sh
RUSTFLAGS=-g cargo build --release --benches
BENCH_FILE=$(find ./target/release/deps -name "benchmarks-*" | grep -v "*.d*" | tail -n1)
echo "profiling ${BENCH_FILE}"
valgrind --tool=cachegrind --I1=32768,8,64 --LL=8388608,16,64 --cachegrind-out-file=cachegrind.out "${BENCH_FILE}" --bench
echo "[done]"