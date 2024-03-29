#include <iterator>
#include "term.h"

namespace groebner_basis {

template <typename Field>
struct NotConstructedTerm {

    Field coefficient;
    std::initializer_list<Degree> degrees;
};

template <typename Field>
static std::vector<Term<Field>> ReduceSimilarInVectorForPolynom(std::vector<Term<Field>> data) {

    for (auto it = data.begin(); it != data.end(); ++it) {
        while ((it + 1) != data.end() && IsDeegreesEqual(*it, *(it + 1))) {
            *it = Term<Field>(it->GetCoefficient() + (it + 1)->GetCoefficient(), *it);
            data.erase(it + 1);
        }
        if (it->GetCoefficient() == 0) {
            it--;
            data.erase(it + 1);
        }
    }
    return data;
}

template <typename Field>
static std::vector<Term<Field>> CreateVectorFromInitializerList(
    const std::initializer_list<NotConstructedTerm<Field>>& terms_list) {

    std::vector<Term<Field>> data;
    data.reserve(terms_list.size());

    for (auto it = terms_list.begin(); it != terms_list.end(); ++it) {
        data.emplace_back(it->coefficient, it->degrees);
    }

    return data;
}

template <typename Field, typename Order>
static std::vector<Term<Field>> CreateOrderedVector(
    const std::initializer_list<NotConstructedTerm<Field>>& terms_list) {

    std::vector<Term<Field>> data = CreateVectorFromInitializerList(terms_list);
    std::sort(data.begin(), data.end(), Order());
    return ReduceSimilarInVectorForPolynom<Field>(data);
}

template <typename Field, typename Order = LexOrder>
class Polynom {
public:
    Polynom(const std::initializer_list<NotConstructedTerm<Field>>& terms_list)
        : data_(CreateOrderedVector<Field, Order>(terms_list)) {
    }

    Polynom(const Term<Field>& term) : data_(ReduceSimilarInVectorForPolynom<Field>({term})) {
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

    size_t TermsCount() const {
        return data_.size();
    }

    bool IsZero() const {
        return data_.empty();
    }

    Polynom operator-() const {
        std::vector<Term<Field>> data;
        data.reserve(data_.size());

        for (const auto& t : (*this)) {
            data.emplace_back(-t);
        }

        return Polynom(data);
    }

    Polynom operator*(const Polynom& other) const {
        std::vector<Term<Field>> result;
        for (const auto& t1 : (*this)) {
            for (const auto& t2 : other) {
                result.push_back(t1 * t2);
            }
        }
        std::sort(result.begin(), result.end(), Order());
        return Polynom(ReduceSimilarInVectorForPolynom<Field>(result));
    }

    Polynom operator+(const Polynom& other) const {
        std::vector<Term<Field>> result;

        std::merge(this->begin(), this->end(), other.begin(), other.end(),
                   std::back_inserter(result), Order());
        return Polynom(ReduceSimilarInVectorForPolynom<Field>(result));
    }

    Polynom operator-(const Polynom& other) const {
        return (*this) + (-other);
    }

private:
    Polynom(const std::vector<Term<Field>>& prepared_vec) : data_(prepared_vec) {
    }

    const std::vector<Term<Field>> data_;
};

template <typename Field, typename Order>
Polynom<Field, Order> operator*(Term<Field> t, const Polynom<Field, Order>& p) {
    return p * t;
}

template <typename Field, typename Order>
Polynom<Field, Order> operator+(Term<Field> t, const Polynom<Field, Order>& p) {
    return p + t;
}

template <typename Field, typename Order>
Polynom<Field, Order> operator-(Term<Field> t, const Polynom<Field, Order>& p) {
    return -(p - t);
}

template <typename Stream, typename Field, typename Order>
Stream& operator<<(Stream& stream, const Polynom<Field, Order>& poly) {
    if (poly.IsZero()) {
        stream << "0";
    } else {
        stream << poly.GetFirstTerm();
        for (auto it = poly.begin() + 1; it != poly.end(); ++it) {
            stream << " + " << *it;
        }
    }
    return stream;
}

}  // namespace groebner_basis
