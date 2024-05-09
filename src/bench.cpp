#include "groebner_basis.h"
#include <algorithm>
#include <boost/rational.hpp>
#include <iostream>
#include <vector>
#include <benchmark/benchmark.h>

namespace bm = benchmark;
namespace gb = groebner_basis;

using Fraction = boost::rational<int64_t>;

static void Cyclic(bm::State &state) {

    size_t n = state.range(0);

    gb::PolynomialsSet<Fraction> s;

    for (size_t i = 1; i < n; ++i) {
        // gb::Polynom<Fraction>::Builder builder;

        std::vector<gb::Monom::Degree> degrees(n, 0);
        std::fill(degrees.begin(), degrees.begin() + i, 1);

        gb::Polynom<Fraction> poly;

        for (size_t j = 0; j < n; ++j) {

            poly = poly + gb::Term<Fraction>(1, gb::BuildMonomFromVectorDegrees(degrees));
            degrees[j] = 0;
            degrees[(j + i) % n] = 1;
        }

        s.Add(poly);
    }

    std::vector<gb::Monom::Degree> degrees(n, 1);
    gb::Polynom<Fraction> poly = gb::Term<Fraction>(1, gb::BuildMonomFromVectorDegrees(degrees));
    poly = poly + gb::Term<Fraction>(-1, {});

    s.Add(poly);

    for (auto _ : state) {
        auto temp = s;
        temp.BuildGreobnerBasis();
        bm::DoNotOptimize(temp);
        std::cout << "begin\n";
        for (auto &p : temp) {
            std::cout << p << "\n";
        }
        std::cout << "end\n";
    }
}

BENCHMARK(Cyclic)->Arg(4);

// BENCHMARK(Cyclic)->Arg(7);

// BENCHMARK(Cyclic)->Arg(8);

BENCHMARK_MAIN();
