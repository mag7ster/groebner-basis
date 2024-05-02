#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include "term.h"

namespace groebner_basis {

template <typename Field, typename Order = LexOrder>
class Polynom {
public:
    class Builder {
    public:
        Builder&& AddTerm(const Field& coef,
                          const std::initializer_list<Monom::Degree>& degrees_list) {

            raw_data_.emplace_back(coef, degrees_list);
            return std::move(*this);
        }

        friend class Polynom;

    private:
        std::vector<Term<Field>>&& GetMovedRawData() {
            return std::move(raw_data_);
        }

        std::vector<Term<Field>> raw_data_;
    };

    Polynom() {  // почему без этого не компилируется???
    }

    Polynom(Builder&& builder)
        : data_(std::make_shared<const std::vector<Term<Field>>>(
              OrderAndReduceVector(std::move(builder.GetMovedRawData())))) {
    }

    Polynom(const Term<Field>& term)
        : data_(std::make_shared<const std::vector<Term<Field>>>(ReduceSimilar({term}))) {
    }

    const Term<Field>& GetLargestTerm() const {
        return data_->front();
    }

    auto begin() const {  // NOLINT
        return data_->begin();
    }

    auto end() const {  // NOLINT
        return data_->end();
    }

    size_t TermsCount() const {
        return data_->size();
    }

    bool IsZero() const {
        return data_->empty();
    }

    Polynom operator-() const {
        std::vector<Term<Field>> data;
        data.reserve(data_->size());

        for (const auto& t : (*this)) {
            data.emplace_back(-t);
        }

        return Polynom(std::move(data));
    }

    friend Polynom operator*(const Polynom& first, const Polynom& second) {
        std::vector<Term<Field>> result;
        result.reserve(first.TermsCount() * second.TermsCount());

        for (const auto& t1 : first) {
            for (const auto& t2 : second) {
                result.push_back(t1 * t2);
            }
        }

        return Polynom(OrderAndReduceVector(std::move(result)));
    }

    friend Polynom operator+(const Polynom& first, const Polynom& second) {
        std::vector<Term<Field>> result;
        result.reserve(first.TermsCount() + second.TermsCount());

        std::merge(first.begin(), first.end(), second.begin(), second.end(),
                   std::back_inserter(result), Order());
        return Polynom(ReduceSimilar(std::move(result)));
    }

    friend Polynom operator-(const Polynom& first, const Polynom& second) {
        return first + (-second);
    }

    bool operator==(const Polynom& other) {
        return *data_ == *other.data_;
    }

    bool operator!=(const Polynom& other) {
        return !(*this == other);
    }

    std::optional<Polynom<Field, Order>> ElementaryReduceBy(const Polynom<Field, Order>& g) const {

        std::optional<Term<Field>> res = FindDivisibleTerm(g.GetLargestTerm());
        if (!res) {
            return std::nullopt;
        }

        Term<Field> divisible_term = res.value();
        Term<Field> t = divisible_term / g.GetLargestTerm();

        return *this - t * g;
    }

    std::optional<Polynom<Field, Order>> ElementaryReduceWithRepeatBy(
        const Polynom<Field, Order>& g) const {

        std::optional<Polynom<Field, Order>> temp = this->ElementaryReduceBy(g);
        if (!temp) {
            return std::nullopt;
        }

        Polynom<Field, Order> res;
        while (temp) {
            res = temp.value();
            temp = res.ElementaryReduceBy(g);
        }

        return res;
    }

private:
    static std::vector<Term<Field>> ReduceSimilar(std::vector<Term<Field>>&& data) {

        // std::unique в итоге не успользуется, потому что он портит элементы

        auto behind_unique = data.begin();
        for (auto current = data.begin() + 1; current != data.end(); ++current) {
            if (current->GetMonom() == behind_unique->GetMonom()) {
                if (current != behind_unique) {
                    *behind_unique =
                        Term(behind_unique->GetCoefficient() + current->GetCoefficient(),
                             behind_unique->GetMonom());
                }
            } else {
                std::swap(*(++behind_unique), *current);
            }
        }
        ++behind_unique;

        data.erase(behind_unique, data.end());

        data.erase(
            std::remove_if(data.begin(), data.end(),
                           [](const Term<Field>& term) { return term.GetCoefficient() == 0; }),
            data.end());

        return std::move(data);
    }

    static std::vector<Term<Field>> OrderAndReduceVector(std::vector<Term<Field>>&& data) {

        std::sort(data.begin(), data.end(), Order());
        return ReduceSimilar(std::move(data));
    }

    Polynom(std::vector<Term<Field>>&& prepared_vec)
        : data_(std::make_shared<const std::vector<Term<Field>>>(std::move(prepared_vec))) {
    }

    std::optional<Term<Field>> FindDivisibleTerm(const Term<Field>& divisor) const {
        for (const auto& t : (*this)) {
            if (t.IsDivisibleBy(divisor)) {
                return t;
            }
        }
        return std::nullopt;
    }

    std::shared_ptr<const std::vector<Term<Field>>> data_;
};

template <typename Stream, typename Field, typename Order>
Stream& operator<<(Stream& stream, const Polynom<Field, Order>& poly) {
    if (poly.IsZero()) {
        stream << "0";
    } else {
        stream << poly.GetLargestTerm();
        for (auto it = poly.begin() + 1; it != poly.end(); ++it) {
            stream << " + " << *it;
        }
    }
    return stream;
}

}  // namespace groebner_basis
