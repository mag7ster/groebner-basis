#include "polynom.h"

namespace groebner_basis {

template <typename Field, typename Order>
Polynom<Field, Order> SPolynom(const Polynom<Field, Order>& f1, const Polynom<Field, Order>& f2) {

    auto lcm = LCM(f1.GetLargestTerm(), f2.GetLargestTerm());
    Term<Field> t1(f2.GetLargestTerm().GetCoefficient(), lcm / f1.GetLargestTerm().GetMonom()),
        t2(f1.GetLargestTerm().GetCoefficient(), lcm / f2.GetLargestTerm().GetMonom());

    return f1 * t1 - f2 * t2;
}

}  // namespace groebner_basis
