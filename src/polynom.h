#include <algorithm>
#include <cstddef>
#include <iterator>
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

        std::vector<Term<Field>>&& GetMovedRawData() {
            return std::move(raw_data_);
        }

    private:
        std::vector<Term<Field>> raw_data_;
    };

    Polynom(Builder&& builder) : data_(OrderAndReduceVector(std::move(builder.GetMovedRawData()))) {
    }

    Polynom(const Term<Field>& term) : data_(ReduceSimilar({term})) {
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

        return Polynom(std::move(data));
    }

    friend Polynom operator*(const Polynom& first, const Polynom& second) {
        std::vector<Term<Field>> result;
        for (const auto& t1 : first) {
            for (const auto& t2 : second) {
                result.push_back(t1 * t2);
            }
        }
        std::sort(result.begin(), result.end(), Order());
        return Polynom(ReduceSimilar(std::move(result)));
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

private:
    static std::vector<Term<Field>> ReduceSimilar(std::vector<Term<Field>>&& data) {

        data.erase(std::unique(data.begin(), data.end(),
                               [](Term<Field>& first, const Term<Field>& second) {
                                   if (first.GetMonom() == second.GetMonom()) {
                                       first =
                                           Term(first.GetCoefficient() + second.GetCoefficient(),
                                                first.GetMonom());
                                       return true;
                                   }
                                   return false;
                               }),
                   data.end());  // Я не уверен, что это не уб. на cppreference сказано только что
                                 // если предикат не устанавливает отношения эквивалентности то
                                 // поведение неопределено. у меня отношение эквивалентности
                                 // сохраняется, хоть и коэффициенты меняются

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

    Polynom(std::vector<Term<Field>>&& prepared_vec) : data_(std::move(prepared_vec)) {
    }

    const std::vector<Term<Field>> data_;
};

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
