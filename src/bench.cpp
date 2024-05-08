#include "groebner_basis.h"
#include <boost/rational.hpp>
#include <cstddef>

#include <benchmark/benchmark.h>

namespace bm = benchmark;
namespace gb = groebner_basis;

using Fraction = boost::rational<int64_t>;

static void Cyclic(bm::State &state) {

    size_t n = state.range(0);

    gb::PolynomialsSet<Fraction> s;

    for (size_t i = 1; i < n; ++i) {
        gb::Polynom<Fraction>::Builder builder;
    }

    for (auto _ : state) {
    }
}

BENCHMARK(Cyclic)->Arg(7);

BENCHMARK(Cyclic)->Arg(8);

BENCHMARK_MAIN();
