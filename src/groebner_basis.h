#include "polynomialsset.h"
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

    std::optional<Term<Field>> res = FindDivisibleTerm(f, g.GetLargestTerm());
    if (!res) {
        return std::nullopt;
    }

    Term<Field> divisible_term = res.value();
    Term<Field> t = divisible_term / g.GetLargestTerm();

    return f - t * g;
}

template <typename Field, typename Order>
Polynom<Field, Order> ReductionByPolynomialsSet(Polynom<Field, Order> f,
                                                const PolynomialsSet<Field, Order> &set) {

    size_t success_reductions;

    do {

        success_reductions = 0;

        for (const auto &g : set) {
            auto temp = ElementaryReduction(f, g);
            while (temp) {
                ++success_reductions;
                f = temp.value();
                temp = ElementaryReduction(f, g);
            }
        }

    } while (success_reductions);

    return f;
}

}  // namespace groebner_basis
