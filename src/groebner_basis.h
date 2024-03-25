#include "polynom.h"

namespace groebner_basis {

template <typename Field, typename Order>
static std::optional<Term<Field>> FindDivisibleTerm(const Polynom<Field, Order> &f,
                                                    const Term<Field> &divisor) {
    for (const auto &t : f.data_) {
        if (t.IsDivisibleBy(divisor)) {
            return t;
        }
    }
    return std::nullopt;
}

template <typename Field, typename Order>
std::optional<Polynom<Field, Order>> ElementaryReduction(Polynom<Field, Order> f,
                                                         Polynom<Field, Order> g) {

    std::optional<Term<Field>> res = FindDivisibleTerm<Field>(f, g.GetFirstTerm()).value();
    if (res == std::nullopt) {
        return std::nullopt;
    }
    Term<Field> divisible_term = res.value();

    //
    // NOT COMPLETED
    //

    return std::nullopt;
}

}  // namespace groebner_basis
