#include "groebner_basis.h"
#include <boost/rational.hpp>
#include <boost/safe_numerics/safe_integer.hpp>
#include <boost/safe_numerics/safe_integer_range.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <cstddef>
#include <iostream>

using namespace boost::safe_numerics;

// using MyInt = boost::safe_numerics::safe_signed_range<-100, 100>;
using MyInt = safe<std::int32_t>;
using Fraction = boost::rational<MyInt>;

namespace gb = groebner_basis;

void PrintCyclic(size_t n) {

    gb::PolynomialsSet<Fraction, gb::GrevLexOrder> s;

    for (size_t i = 1; i < n; ++i) {
        // gb::Polynom<Fraction>::Builder builder;

        std::vector<gb::Monom::Degree> degrees(n, 0);
        std::fill(degrees.begin(), degrees.begin() + i, 1);

        gb::Polynom<Fraction, gb::GrevLexOrder> poly;

        for (size_t j = 0; j < n; ++j) {

            poly = poly + gb::Term<Fraction>(Fraction(1), gb::BuildMonomFromVectorDegrees(degrees));
            degrees[j] = 0;
            degrees[(j + i) % n] = 1;
        }

        s.Add(poly);
    }

    std::vector<gb::Monom::Degree> degrees(n, 1);
    gb::Polynom<Fraction, gb::GrevLexOrder> poly =
        gb::Term<Fraction>(Fraction(1), gb::BuildMonomFromVectorDegrees(degrees));
    poly = poly + gb::Term<Fraction>(Fraction(-1), {});

    s.Add(poly);

    std::cout << "Cyclic " << n << " \n";
    for (auto &p : s) {
        std::cout << p << "\n";
    }
    std::cout << "\n";

    try {

        s.BuildGreobnerBasis();
        for (auto &p : s) {
            std::cout << p << "\n";
        }
        std::cout << "\n";

    } catch (...) {
        std::cerr << boost::current_exception_diagnostic_information();
    }
}

int main() {
    // groebner_basis::Polynom<Fraction> a =
    //     groebner_basis::Polynom<Fraction>::Builder().AddTerm(1, {2, 3}).AddTerm(4, {5,
    //     6}).AddTerm(
    //         4, {5, 6});
    // groebner_basis::Monom m1({1, 2, 3, 0, 0});
    // groebner_basis::Term<Fraction> t1(1, {2, 3}), t2(2, {3, 4});
    // groebner_basis::Polynom<Fraction> b = t1;
    // std::cout << m1 << "\n";
    // std::cout << t1 << "\n";
    // std::cout << t2 << "\n";
    // std::cout << groebner_basis::LexOrder()(t1, t2) << "\n";
    // std ::cout << t2.IsDivisibleBy(t1) << "\n";

    // std ::cout << a << "\n";
    // std ::cout << b << "\n";
    // std::cout << a + b << "\n";

    // std::cout << a * b << "\n";
    // std::cout << a - b << "\n";

    // std::cout << a * t2 << "\n";
    // std::cout << t2 * a << "\n";

    // std::cout << a + t2 << "\n";
    // std::cout << t2 + a << "\n";

    // std::cout << a - t2 << "\n";
    // std::cout << t2 - a << "\n";

    // std::cout << "Self -\n";
    // std::cout << a - a << "\n";

    // std::cout << "Elem\n";
    // auto res = a.ElementaryReduceBy(a);
    // std::cout << res.value() << "\n";

    // groebner_basis::PolynomialsSet g({a, b});
    // auto res1 = groebner_basis::PolynomialsSet({a, b}).Reduce(a);
    // std::cout << res1.value() << "\n";

    // std::cout << t1.GetMonom() << " | " << t2.GetMonom() << "\n";
    // auto lol = groebner_basis::LCM(t1.GetMonom(), t2.GetMonom());
    // std::cout << lol << "\n" << lol / t1.GetMonom() << "\n" << lol / t2.GetMonom() << "\n";

    // std::cout << a << " | " << b << "\n";
    // std::cout << groebner_basis::SPolynom(a, b) << "\n";

    // g.BuildGreobnerBasis();
    // std::cout << "Greobner\n";
    // for (const auto& f : g) {
    //     std::cout << f << "\n";
    // }

    PrintCyclic(4);
    PrintCyclic(5);

    // for (int i = 2; i < 7; ++i) {
    //     PrintCyclic(i);
    // }

    return 0;
}
