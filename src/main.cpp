#include "groebner_basis.h"
#include <boost/rational.hpp>
#include <iostream>

using Fraction = boost::rational<int64_t>;

int main() {
    groebner_basis::Polynom<Fraction> a =
        groebner_basis::Polynom<Fraction>::Builder().AddTerm(1, {2, 3}).AddTerm(4, {5, 6}).AddTerm(
            4, {5, 6});
    groebner_basis::Monom m1({1, 2, 3, 0, 0});
    groebner_basis::Term<Fraction> t1(1, {2, 3}), t2(2, {3, 4});
    groebner_basis::Polynom<Fraction> b = t1;
    std::cout << m1 << "\n";
    std::cout << t1 << "\n";
    std::cout << t2 << "\n";
    std::cout << groebner_basis::LexOrder()(t1, t2) << "\n";
    std ::cout << t2.IsDivisibleBy(t1) << "\n";

    std ::cout << a << "\n";
    std ::cout << b << "\n";
    std::cout << a + b << "\n";

    std::cout << a * b << "\n";
    std::cout << a - b << "\n";

    std::cout << a * t2 << "\n";
    std::cout << t2 * a << "\n";

    std::cout << a + t2 << "\n";
    std::cout << t2 + a << "\n";

    std::cout << a - t2 << "\n";
    std::cout << t2 - a << "\n";

    std::cout << "Self -\n";
    std::cout << a - a << "\n";

    std::cout << "Elem\n";
    auto res = a.ElementaryReduceBy(a);
    std::cout << res.value() << "\n";

    groebner_basis::PolynomialsSet g({a, b});
    auto res1 = groebner_basis::PolynomialsSet({a, b}).Reduce(a);
    std::cout << res1.value() << "\n";
}
