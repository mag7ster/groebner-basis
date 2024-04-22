#include "polynom.h"

namespace groebner_basis {

template <typename Field, typename Order>
Polynom<Field, Order> SPolynom(const Polynom<Field, Order>& f, const Polynom<Field, Order>& g) {
    return Polynom<Field, Order>();
}

}  // namespace groebner_basis
