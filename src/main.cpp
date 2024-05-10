#include "groebner_basis.h"
#include <boost/rational.hpp>
#include <boost/safe_numerics/safe_integer.hpp>
#include <cassert>
#include <fstream>
#include <iostream>

// using MyInt = boost::safe_numerics::safe_signed_range<-100, 100>;
using MyInt = boost::safe_numerics::safe<std::int64_t>;
// using MyInt = std::int32_t;
using Fraction = boost::rational<MyInt>;

namespace gb = groebner_basis;

// void PrintCyclic(size_t n) {

//     gb::PolynomialsSet<Fraction, gb::GrevLexOrder> s;

//     for (size_t i = 1; i < n; ++i) {
//         // gb::Polynom<Fraction>::Builder builder;

//         std::vector<gb::Monom::Degree> degrees(n, 0);
//         std::fill(degrees.begin(), degrees.begin() + i, 1);

//         gb::Polynom<Fraction, gb::GrevLexOrder> poly;

//         for (size_t j = 0; j < n; ++j) {

//             poly = poly + gb::Term<Fraction>(Fraction(1),
//             gb::BuildMonomFromVectorDegrees(degrees)); degrees[j] = 0; degrees[(j + i) % n] = 1;
//         }

//         s.Add(poly);
//     }

//     std::vector<gb::Monom::Degree> degrees(n, 1);
//     gb::Polynom<Fraction, gb::GrevLexOrder> poly =
//         gb::Term<Fraction>(Fraction(1), gb::BuildMonomFromVectorDegrees(degrees));
//     poly = poly + gb::Term<Fraction>(Fraction(-1), {});

//     s.Add(poly);

//     std::cout << "Cyclic " << n << " \n";
//     for (auto& p : s) {
//         std::cout << p << "\n";
//     }
//     std::cout << "\n";

//     try {

//         s.BuildGreobnerBasis();
//         for (auto& p : s) {
//             std::cout << p << "\n";
//         }
//         std::cout << "\n";

//     } catch (...) {
//         std::cerr << boost::current_exception_diagnostic_information();
//     }
// }

void Check(gb::PolynomialsSet<Fraction, gb::GrevLexOrder>& find,
           gb::PolynomialsSet<Fraction, gb::GrevLexOrder>& ans) {

    std::cout << "test\n";
    for (auto& f : find) {
        std::cout << f << "\n";
    }

    find.BuildGreobnerBasis();
    std::sort(ans.begin(), ans.end());

    std::cout << "test groebner\n";
    for (auto& f : find) {
        std::cout << f << "\n";
    }

    std::cout << "ans\n";
    for (auto& f : ans) {
        std::cout << f << "\n";
    }
    assert(find == ans);
}

int main() {

    int state;
    std::string str;

    gb::PolynomialsSet<Fraction, gb::GrevLexOrder> find, ans;

    std::ifstream file("/home/mukar/groebner-basis/tests.txt");

    while (file >> str) {
        // std::cout << str << "\n";

        if (str == "test") {
            state = 0;
            continue;
        }
        if (str == "Groebner") {
            state = 1;
            continue;
        }
        if (str == "calc") {
            Check(find, ans);
            find.Clear();
            ans.Clear();
            std::cout << "A\n";
            continue;
        }

        std::stringstream ss;
        ss << str;

        gb::Polynom<Fraction, gb::GrevLexOrder> poly;

        Fraction frac = 1;
        std::vector<gb::Monom::Degree> degs(3, 0);
        int sign = 1;

        if (ss.peek() == '-') {
            frac = 0;
        }

        while (ss.good()) {
            if (ss.peek() == '-') {
                poly =
                    poly + gb::Term<Fraction>(sign * frac, gb::BuildMonomFromVectorDegrees(degs));

                ss.get();

                frac = 1;
                degs = {0, 0, 0};
                sign = -1;

            } else if (ss.peek() == '+') {
                poly =
                    poly + gb::Term<Fraction>(sign * frac, gb::BuildMonomFromVectorDegrees(degs));

                ss.get();

                frac = 1;
                degs = {0, 0, 0};
                sign = 1;

            } else if (ss.peek() == 'x' || ss.peek() == 'y' || ss.peek() == 'z') {

                char c = ss.get();

                gb::Monom::Degree deg = 1;

                if (ss.peek() == '^') {
                    ss.get();
                    deg = ss.get();
                    deg -= '0';
                }

                degs[c - 'x'] = deg;

            } else if (std::isdigit(ss.peek())) {
                ss >> frac;
            }
        }

        poly = poly + gb::Term<Fraction>(sign * frac, gb::BuildMonomFromVectorDegrees(degs));

        if (state == 0) {
            find.Add(poly);
        } else {
            ans.Add(poly);
        }
    }
    file.close();

    // PrintCyclic(4);
    // PrintCyclic(5);

    return 0;
}
