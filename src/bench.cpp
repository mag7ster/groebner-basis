#include "groebner_basis.h"
#include <algorithm>
#include <cstdint>
#include <vector>
#include <benchmark/benchmark.h>
#include "types.h"

namespace bm = benchmark;
namespace gb = groebner_basis;

using ModInt = gb::Modulus<std::int64_t, 998244353>;

static gb::PolynomialsSet<ModInt> BuildCyclic(int n) {
    gb::PolynomialsSet<ModInt> s;

    for (size_t i = 1; i < n; ++i) {

        std::vector<gb::Monom::Degree> degrees(n, 0);
        std::fill(degrees.begin(), degrees.begin() + i, 1);

        gb::Polynom<ModInt>::Builder poly;

        for (size_t j = 0; j < n; ++j) {

            poly = poly.AddTerm(1, gb::Monom::BuildFromVectorDegrees(degrees));
            degrees[j] = 0;
            degrees[(j + i) % n] = 1;
        }

        s.Add(poly.BuildPolynom());
    }

    std::vector<gb::Monom::Degree> degrees(n, 1);
    gb::Polynom<ModInt>::Builder poly;
    poly = poly.AddTerm(1, gb::Monom::BuildFromVectorDegrees(degrees));
    poly = poly.AddTerm(-1, {});

    s.Add(poly.BuildPolynom());
    return s;
}

static void Cyclic(bm::State &state) {

    size_t n = state.range(0);

    auto s = BuildCyclic(n);

    for (auto _ : state) {
        auto temp = s;
        temp.BuildGreobnerBasis();
        bm::DoNotOptimize(temp);
    }
}

BENCHMARK(Cyclic)->Arg(4)->Iterations(1000)->Unit(bm::kMillisecond);
BENCHMARK(Cyclic)->Arg(5)->Iterations(10)->Unit(bm::kMillisecond);
BENCHMARK(Cyclic)->Arg(6)->Iterations(1)->Unit(bm::kSecond);

BENCHMARK_MAIN();
