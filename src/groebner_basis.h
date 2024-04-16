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
std::optional<Polynom<Field, Order>> ElementaryReductionWithRepeat(const Polynom<Field, Order> &f,
                                                                   const Polynom<Field, Order> &g) {

    std::optional<Polynom<Field, Order>> temp = ElementaryReduction(f, g);
    if (!temp) {
        return std::nullopt;
    }

    Polynom<Field, Order> res;
    while (temp) {
        res = temp.value();
        temp = ElementaryReduction(res, g);
    }

    return res;
}

template <typename Field, typename Order>
static std::optional<Polynom<Field, Order>> TryReductionForOnePass(
    const Polynom<Field, Order> &f, const PolynomialsSet<Field, Order> &set) {

    Polynom<Field, Order> res = f;
    bool success = false;

    for (auto &g : set) {
        auto temp = ElementaryReductionWithRepeat(res, g);
        if (temp) {
            success = true;
            res = temp.value();
        }
    }

    if (!success) {
        return std::nullopt;
    }

    return res;
}

template <typename Field, typename Order>
std::optional<Polynom<Field, Order>> ReductionByPolynomialsSet(
    const Polynom<Field, Order> &f, const PolynomialsSet<Field, Order> &set) {

    std::optional<Polynom<Field, Order>> temp = TryReductionForOnePass(f, set);
    if (!temp) {
        return std::nullopt;
    }

    Polynom<Field, Order> res;
    while (temp) {
        res = temp.value();
        temp = TryReductionForOnePass(res, set);
    }

    return res;
}

// not completed

template <typename Field, typename Order>
bool TryAutoReductionForOnePass(PolynomialsSet<Field, Order> *set) {
    for (auto it = set->begin(); it != set->end(); ++it) {
        auto f = *it;
        set->Erase(it);
        auto temp = ReductionByPolynomialsSet(f, set);
        if (!temp) {
            set->AddAt(it, f);
        } else if (temp.value() != Term<Field>(0)) {
            set->AddAt(it, f);
        }
    }
}

template <typename Field, typename Order>
void AutoReduction(PolynomialsSet<Field, Order> *set) {
}

}  // namespace groebner_basis
