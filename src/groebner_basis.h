#include <cstddef>
#include <vector>
#include "functions.h"

namespace groebner_basis {

template <typename Field, typename Order = GrevLexOrder>
class PolynomialsSet {

    using Polynom = Polynom<Field, Order>;

    using Container = std::vector<Polynom>;
    using Iterator = typename Container::iterator;
    using ConstIterator = typename Container::const_iterator;

public:
    PolynomialsSet(std::initializer_list<Polynom> poly_list) : data_(poly_list) {
    }

    PolynomialsSet() = default;

    Iterator begin() {  // NOLINT
        return data_.begin();
    }

    ConstIterator begin() const {  // NOLINT
        return data_.begin();
    }

    Iterator end() {  // NOLINT
        return data_.end();
    }

    ConstIterator end() const {  // NOLINT
        return data_.end();
    }

    size_t Size() const {
        return data_.size();
    }

    bool operator==(const PolynomialsSet &other) const {
        return data_ == other.data_;
    }

    bool operator!=(const PolynomialsSet &other) const {
        return !(*this == other);
    }

    void Add(const Polynom &poly) {
        if (poly.IsZero()) {
            return;
        }

        data_.emplace_back(poly);
        data_.back() =
            data_.back() * Term<Field>(Field(1) / data_.back().GetLargestTerm().GetCoefficient());
    }

    void Add(Polynom &&poly) {
        if (poly.IsZero()) {
            return;
        }

        data_.emplace_back(std::move(poly));
        data_.back() =
            data_.back() * Term<Field>(Field(1) / data_.back().GetLargestTerm().GetCoefficient());
    }

    void Erase(Iterator it) {
        if (it != data_.end()) {
            std::swap(*it, data_.back());
            data_.pop_back();
        }
    }

    void Erase(const Polynom &poly) {
        auto it = std::find(data_.begin(), data_.end(), poly);
        Erase(it);
    }

    void Clear() {
        data_.clear();
    }

    std::optional<Polynom> Reduce(const Polynom &f) const {

        std::optional<Polynom> temp = TryReductionForOnePass(f);
        if (!temp) {
            return std::nullopt;
        }

        Polynom res;
        while (temp) {
            res = temp.value();
            temp = TryReductionForOnePass(res);
        }

        return res;
    }

    void AutoReduction() {
        auto it = begin();
        for (size_t i = 0; i < Size(); ++i, ++it) {

            auto f = *it;
            this->Erase(it);
            auto temp = this->Reduce(f);
            if (!temp) {
                this->AddAt(it, f);
            } else {
                if (!temp.value().IsZero()) {
                    this->AddAt(it, temp.value());
                } else {
                    it--;
                    i--;
                }
            }
        }

        for (auto &f : (*this)) {
            f = f * Term<Field>(Field(1) / f.GetLargestTerm().GetCoefficient());
        }
    }

    void BuildGreobnerBasis() {

        BuildUnReducedGroebnerBasis();
        Minimize();
        AutoReduction();
        std::sort(begin(), end());
    }

private:
    void AddAt(Iterator it, const Polynom &poly) {
        Add(poly);
        std::swap(*it, data_.back());
    }

    void Minimize() {

        for (auto it1 = begin(); it1 != end(); ++it1) {
            bool remove = false;
            std::swap(*it1, *begin());

            for (auto it2 = begin() + 1; it2 != end(); ++it2) {
                if (data_.front().GetLargestTerm().IsDivisibleBy(it2->GetLargestTerm())) {
                    remove = true;
                }
            }

            if (remove) {
                Erase(begin());
                it1--;
            }
        }
    }

    std::optional<Polynom> TryReductionForOnePass(const Polynom &f) const {

        Polynom res = f;
        bool success = false;

        for (auto &g : (*this)) {

            auto temp = res.ElementaryReduceWithRepeatBy(g);
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

    void BuildUnReducedGroebnerBasis() {

        for (size_t i = 0; i < Size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                auto s = SPolynom(data_[i], data_[j]);

                if (s.IsZero()) {
                    continue;
                }

                auto r_ij = Reduce(s);

                if (!r_ij) {
                    Add(s);
                } else if (!r_ij.value().IsZero()) {
                    Add(r_ij.value());
                }
            }
        }
    }

    Container data_;
};

}  // namespace groebner_basis
