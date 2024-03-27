
#include "term.h"

namespace groebner_basis {

template <typename Field>
struct NotConstructedTerm {

    Field coefficient;
    std::initializer_list<Degree> degrees;
};

template <typename Field, typename Order>
static std::vector<Term<Field>> CreatePolynom(
    const std::initializer_list<NotConstructedTerm<Field>>& terms_list) {

    std::vector<Term<Field>> data;
    data.reserve(terms_list.size());

    for (auto it = terms_list.begin(); it != terms_list.end(); ++it) {
        data.emplace_back(it->coefficient, it->degrees);
    }

    std::sort(data.begin(), data.end(), Order());

    for (auto it = data.begin(); it != data.end(); ++it) {
        for (auto next = it + 1; next != data.end() && IsDeegreesEqual(*it, *next); ++next) {
            *it = Term<Field>(it->GetCoefficient() + next->GetCoefficient(), *it);
            data.erase(next);
        }
    }

    return data;
}

template <typename Field, typename Order = LexOrder>
class Polynom {
public:
    Polynom(const std::initializer_list<NotConstructedTerm<Field>>& terms_list)
        : data_(CreatePolynom<Field, Order>(terms_list)) {
    }

    const Term<Field>& GetFirstTerm() const {
        return data_.front();
    }

    auto begin() const {  // NOLINT
        return data_.begin();
    }

    auto end() const {  // NOLINT
        return data_.end();
    }

private:
    std::vector<Term<Field>> data_;
};

template <typename Stream, typename Field, typename Order>
Stream& operator<<(Stream& stream, const Polynom<Field, Order>& poly) {
    stream << poly.GetFirstTerm();
    for (auto it = poly.begin() + 1; it != poly.end(); ++it) {
        stream << " + " << *it;
    }
    return stream;
}

}  // namespace groebner_basis
