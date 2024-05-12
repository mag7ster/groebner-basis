#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <optional>
#include "term.h"
#include <sstream>

namespace groebner_basis {

template <typename Field, typename Order = GrevLexOrder>
class Polynom {
public:
    using Term = Term<Field>;

    class Builder {
    public:
        Builder& AddTerm(const Field& coef, std::initializer_list<Monom::Degree> degrees_list) {

            raw_data_.emplace_back(coef, degrees_list);
            return *this;
        }

        Builder& AddTerm(const Field& coef, const Monom& monom) {

            raw_data_.emplace_back(coef, monom);
            return *this;
        }

        Polynom BuildPolynom() {
            return Polynom(OrderAndReduceVector(std::move(raw_data_)));
        }

    private:
        std::vector<Term> raw_data_;
    };
    friend class Builder;

    Polynom() = default;

    Polynom(Builder&& builder) : Polynom(builder.BuildPolynom()) {
        assert(IsCorrect());
    }

    Polynom(const Term& term) : Polynom(ReduceSimilar({term})) {
        assert(IsCorrect());
    }

    static Polynom BuildFromString(const std::string& str) {
        return ParseAndBuild(str);
    }

    const Term& GetLargestTerm() const {
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
        std::vector<Term> data;
        data.reserve(data_->size());

        for (const auto& t : (*this)) {
            data.emplace_back(-t);
        }

        return Polynom(std::move(data));
    }

    friend Polynom operator*(const Polynom& first, const Polynom& second) {
        std::vector<Term> result;
        result.reserve(first.TermsCount() * second.TermsCount());

        for (const auto& t1 : first) {
            for (const auto& t2 : second) {
                result.push_back(t1 * t2);
            }
        }

        return Polynom(OrderAndReduceVector(std::move(result)));
    }

    friend Polynom operator+(const Polynom& first, const Polynom& second) {
        std::vector<Term> result;
        result.reserve(first.TermsCount() + second.TermsCount());

        std::merge(first.begin(), first.end(), second.begin(), second.end(),
                   std::back_inserter(result), Order());
        return Polynom(ReduceSimilar(std::move(result)));
    }

    friend Polynom operator-(const Polynom& first, const Polynom& second) {
        return first + (-second);
    }

    friend bool operator==(const Polynom& first, const Polynom& second) {
        return *first.data_ == *second.data_;
    }

    friend bool operator!=(const Polynom& first, const Polynom& second) {
        return !(first == second);
    }

    friend bool operator<(const Polynom& first, const Polynom& second) {

        for (auto it1 = first.begin(), it2 = second.begin();
             it1 != first.end() && it2 != second.end(); ++it1, ++it2) {

            if (*it1 == *it2) {
                continue;
            } else if (it1->GetMonom() == it2->GetMonom()) {
                return it1->GetCoefficient() > it2->GetCoefficient();
            } else {
                return Order()(*it1, *it2);
            }
        }

        return first.TermsCount() > second.TermsCount();
    }

    template <typename Stream>
    friend Stream& operator<<(Stream& stream, const Polynom& poly) {
        if (poly.IsZero()) {
            stream << "0";
        } else {
            stream << poly.GetLargestTerm();
            for (auto it = poly.begin() + 1; it != poly.end(); ++it) {
                if (it->GetCoefficient() < 0) {
                    stream << " - " << Term(-it->GetCoefficient(), it->GetMonom());
                } else {
                    stream << " + " << *it;
                }
            }
        }
        return stream;
    }

    std::optional<Polynom> ElementaryReduceBy(const Polynom& g) const {

        std::optional<Term> optdiv = FindDivisibleTerm(g.GetLargestTerm());
        if (!optdiv) {
            return std::nullopt;
        }

        Term divisible_term = optdiv.value();
        Term t = divisible_term / g.GetLargestTerm();

        return *this - t * g;
    }

    std::optional<Polynom> ElementaryReduceWithRepeatBy(const Polynom& g) const {

        std::optional<Polynom> temp = this->ElementaryReduceBy(g);
        if (!temp) {
            return std::nullopt;
        }

        Polynom res;
        while (temp) {
            res = temp.value();
            temp = res.ElementaryReduceBy(g);
        }

        return res;
    }

private:
    static std::vector<Term> ReduceSimilar(std::vector<Term>&& data) {

        if (data.empty()) {
            return std::move(data);
        }

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
                           [](const Term& term) { return term.GetCoefficient() == Field(0); }),
            data.end());

        return std::move(data);
    }

    static std::vector<Term> OrderAndReduceVector(std::vector<Term>&& data) {

        std::sort(data.begin(), data.end(), Order());
        return ReduceSimilar(std::move(data));
    }

    Polynom(std::vector<Term>&& prepared_vec)
        : data_(std::make_shared<const std::vector<Term>>(std::move(prepared_vec))) {
        assert(IsCorrect());
    }

    static Polynom ParseAndBuild(const std::string& str) {  // works only with xyz
        std::stringstream ss;
        ss << str;

        Polynom::Builder builder;

        Field current_coef = 1;
        std::vector<Monom::Degree> degs(3, 0);
        int sign = 1;

        if (ss.peek() == '-') {
            current_coef = 0;
        }

        while (ss.good()) {
            if (ss.peek() == '-') {
                builder = builder.AddTerm(sign * current_coef, Monom::BuildFromVectorDegrees(degs));

                ss.get();

                current_coef = 1;
                degs = {0, 0, 0};
                sign = -1;

            } else if (ss.peek() == '+') {
                builder = builder.AddTerm(sign * current_coef, Monom::BuildFromVectorDegrees(degs));

                ss.get();

                current_coef = 1;
                degs = {0, 0, 0};
                sign = 1;

            } else if (ss.peek() == 'x' || ss.peek() == 'y' || ss.peek() == 'z') {

                char c = ss.get();

                Monom::Degree deg = 1;

                if (ss.peek() == '^') {
                    ss.get();
                    ss >> deg;
                }

                degs[c - 'x'] = deg;

            } else if (std::isdigit(ss.peek())) {
                int num;
                ss >> num;
                current_coef = num;
            }
        }

        builder = builder.AddTerm(sign * current_coef, Monom::BuildFromVectorDegrees(degs));
        return Polynom(std::move(builder));
    }

    std::optional<Term> FindDivisibleTerm(const Term& divisor) const {

        auto it = std::find_if(this->begin(), this->end(),
                               [&](const Term& t) { return t.IsDivisibleBy(divisor); });

        if (it == this->end()) {
            return std::nullopt;
        }
        return *it;
    }

    bool Find(const Monom& m) const {
        for (const auto& t : (*this)) {
            if (t.GetMonom() == m) {
                return true;
            }
        }
        return false;
    }

    bool IsCorrect() {
        if (IsZero()) {
            return true;
        }

        for (auto it = begin() + 1; it != end(); ++it) {
            if (!Order()(*(it - 1), *it)) {
                return false;
            }
        }
        for (auto it = begin(); it != end(); ++it) {
            for (auto it2 = it + 1; it2 != end(); ++it2) {
                if (it->GetMonom() == it2->GetMonom()) {
                    return false;
                }
            }
        }

        for (auto it = begin(); it != end(); ++it) {
            if (it->GetCoefficient() == Field(0)) {
                return false;
            }
        }

        return true;
    }

    std::shared_ptr<const std::vector<Term>> data_ = std::make_shared<const std::vector<Term>>();
};

}  // namespace groebner_basis
