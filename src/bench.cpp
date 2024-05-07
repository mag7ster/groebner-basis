#include "groebner_basis.h"
#include <boost/rational.hpp>

#include <benchmark/benchmark.h>

namespace bm = benchmark;

static void Cyclic(bm::State &state) {
}

// The compiler will just optimize everything out.
// After the first run, the value of `c` won't change.
// The benchmark will show 0ns per iteration.
BENCHMARK(Cyclic);

BENCHMARK_MAIN();
