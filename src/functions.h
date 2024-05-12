#include <cassert>
#include "polynom.h"

namespace groebner_basis {

inline Monom LCM(const Monom& m1, const Monom& m2) {

    size_t size =
        std::max(m1.FirstIndexAfterLastNonZeroDegree(), m2.FirstIndexAfterLastNonZeroDegree());
    std::vector<Monom::Degree> lcm(size);

    for (size_t i = 0; i < size; ++i) {
        lcm[i] = std::max(m1.Deg(i), m2.Deg(i));
    }

    return Monom(std::move(lcm));
}

template <typename Field, typename Order>
Polynom<Field, Order> SPolynom(const Polynom<Field, Order>& f1, const Polynom<Field, Order>& f2) {

    auto lcm = LCM(f1.GetLargestTerm(), f2.GetLargestTerm());
    Term<Field> t1(f2.GetLargestTerm().GetCoefficient(), lcm / f1.GetLargestTerm().GetMonom()),
        t2(f1.GetLargestTerm().GetCoefficient(), lcm / f2.GetLargestTerm().GetMonom());

    assert(f1.GetLargestTerm() * t1 == f2.GetLargestTerm() * t2);

    return f1 * t1 - f2 * t2;
}

}  // namespace groebner_basis
