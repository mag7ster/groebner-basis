#include <algorithm>
#include <initializer_list>
#include <optional>
#include <string>
#include <vector>
#include "tostringhelper.h"

namespace groebner_basis {

class Monom;

class LexOrder {
public:
    bool operator()(const Monom &a, const Monom &b);
};

class Monom {
public:
    Monom(std::initializer_list<int> list);

    friend class LexOrder;

protected:
    std::vector<int> degs_;
};

template <typename Field>
class Term : public Monom {
public:
    Term(std::pair<Field, std::initializer_list<int>> term)
        : Monom(term.second), coef_(term.first) {
    }

    std::string ToString() {
        std::string res = tostringhelper::ToString<Field>(coef_) + "*";
        for (size_t i = 0; i < degs_.size(); ++i) {
            res += "x" + std::to_string(i) + "^{" + std::to_string(degs_[i]) + "}";
        }
        return res;
    }

    template <typename T>
    friend bool IsDivisible(const groebner_basis::Term<T> &term,
                            const groebner_basis::Term<T> &divisor);

protected:
    Field coef_;
};

template <typename Field>
class Polynom {
public:
    Polynom(std::initializer_list<std::pair<Field, std::initializer_list<int>>> list) {
        data_.reserve(list.size());
        for (auto it = list.begin(); it != list.end(); ++it) {
            data_.emplace_back(*it);
        }
    }

    std::string ToString() {
        std::string res;
        for (auto &term : data_) {
            res += term.ToString() + "+";
        }
        res.pop_back();
        return res;
    }

    Term<Field> GetFirstTerm() {
        return data_.front();
    }

    template <typename T>
    friend std::optional<Term<T>> FindDivisibleTerm(const Polynom<T> &f, const Term<T> &t);

    template <typename T, typename Functor>
    friend std::optional<Polynom<T>> ElementaryReduction(Polynom<T> f, Polynom<T> g, Functor order);

private:
    std::vector<Term<Field>> data_;
};

template <typename Field>
static bool IsDivisible(const Term<Field> &term, const Term<Field> &divisor) {

    for (size_t i = 0; i < std::min(term.degs_.size(), divisor.degs_.size()); ++i) {
        if (divisor.degs_[i] > term.degs_[i]) {
            return false;
        }
    }
    return true;
}

template <typename Field>
static std::optional<Term<Field>> FindDivisibleTerm(const Polynom<Field> &f,
                                                    const Term<Field> &divisor) {
    for (auto &t : f.data_) {
        if (IsDivisible(t, divisor)) {
            return t;
        }
    }
    return std::nullopt;
}

template <typename Field, typename Functor>
std::optional<Polynom<Field>> ElementaryReduction(Polynom<Field> f, Polynom<Field> g,
                                                  Functor order) {
    std::sort(f.data_.begin(), f.data_.end(), order);
    std::sort(g.data_.begin(), g.data_.end(), order);

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