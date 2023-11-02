#include <benchmark/benchmark.h>

#include <functional>
#include <numbers>
#include <vector>

#include "hdnum.hh"

static void MatrixVector(benchmark::State& state) {
    hdnum::Vector<double> v(state.range(0));
    hdnum::SparseMatrix<double> A(state.range(0), state.range(0));
    A = A.matchingIdentity();
    // fill with 1, 2, 3, 4 ... N
    std::generate(std::begin(v), std::end(v),
                  [n = 0]() mutable { return n++; });
    for (auto _ : state) {
        hdnum::Vector<double> res = A * v;
    }
    state.SetComplexityN(state.range(0));
}
BENCHMARK(MatrixVector)
    ->RangeMultiplier(2)
    ->Range(1 << 6, 1 << 22)
    ->Complexity(benchmark::oN);

BENCHMARK_MAIN();
