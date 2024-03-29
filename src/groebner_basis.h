#include "polynom.h"
#include <optional>

namespace groebner_basis {

template <typename Field, typename Order>
static std::optional<Term<Field>> FindDivisibleTerm(const Polynom<Field, Order> &f,
                                                    const Term<Field> &divisor) {
    for (const auto &t : f) {
        if (t.IsDivisibleBy(divisor)) {
            return t;
        }
    }
    return std::nullopt;
}

template <typename Field, typename Order>
std::optional<Polynom<Field, Order>> ElementaryReduction(const Polynom<Field, Order> &f,
                                                         const Polynom<Field, Order> &g) {

    std::optional<Term<Field>> res = FindDivisibleTerm(f, g.GetFirstTerm());
    if (res == std::nullopt) {
        return std::nullopt;
    }

    Term<Field> divisible_term = res.value();
    Term<Field> t = divisible_term / g.GetFirstTerm();

    return f - t * g;
}

}  // namespace groebner_basis
