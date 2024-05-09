#include <cstddef>
#include <iostream>
#include <vector>
#include "functions.h"

namespace groebner_basis {

template <typename Field, typename Order = GrevLexOrder>
class PolynomialsSet {

    using Polynom = Polynom<Field, Order>;

    using Container = std::vector<Polynom>;
    using Iterator = typename Container::iterator;

public:
    PolynomialsSet(const std::initializer_list<Polynom> &poly_list) : data_(poly_list) {
    }

    PolynomialsSet() {
    }

    auto begin() {  // NOLINT
        return data_.begin();
    }

    auto begin() const {  // NOLINT
        return data_.begin();
    }

    auto end() {  // NOLINT
        return data_.end();
    }

    auto end() const {  // NOLINT
        return data_.end();
    }

    size_t Size() const {
        return data_.size();
    }

    void Add(const Polynom &poly) {
        data_.emplace_back(poly);
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

    unsigned int AutoReduction() {
        unsigned int count_passes = 0;
        while (TryAutoReductionForOnePass()) {
            count_passes += 1;
        }
        return count_passes;
    }

    unsigned int BuildGreobnerBasis() {
        unsigned int count_passes = 0;
        while (TryBuildGreobnerBasisForOnePass()) {
            count_passes += 1;

            std::cout << "CUR BASIS " << count_passes << "---------\n";
            for (const auto &f : (*this)) {
                std::cout << f << "\n\n";
            }
            std::cout << "--------------------------\n";
        }

        for (auto &f : (*this)) {
            f = f * Term<Field>(Field(1) / f.GetLargestTerm().GetCoefficient());
        }

        std::sort(begin(), end());

        return count_passes;
    }

private:
    void AddAt(Iterator it, const Polynom &poly) {
        Add(poly);
        std::swap(*it, data_.back());
    }

    std::optional<Polynom> TryReductionForOnePass(const Polynom &f) const {

        Polynom res = f;
        bool success = false;

        int i = 0;
        for (auto &g : (*this)) {

            auto temp = res.ElementaryReduceBy(g);
            if (temp) {
                success = true;
                res = temp.value();
            }
            ++i;
        }

        if (!success) {
            return std::nullopt;
        }

        return res;
    }

    bool TryAutoReductionForOnePass() {

        bool success = false;

        auto it = this->begin();
        for (size_t i = 0; i < Size(); ++i, ++it) {

            auto f = *it;
            this->Erase(it);
            auto temp = this->Reduce(f);
            if (!temp) {
                this->AddAt(it, f);
            } else {
                success = true;
                if (!temp.value().IsZero()) {
                    this->AddAt(it, temp.value());
                }
            }
        }

        return success;
    }

    bool TryBuildGreobnerBasisForOnePass() {
        std::vector<Polynom> r;
        r.reserve(Size() * (Size() - 1) / 2);
        for (auto it1 = begin(); it1 != end(); ++it1) {
            for (auto it2 = it1 + 1; it2 != end(); ++it2) {
                auto s = SPolynom(*it1, *it2);
                auto r_ij = Reduce(s);

                if (r_ij && !r_ij.value().IsZero()) {
                    r.emplace_back(r_ij.value());
                }
            }
        }

        for (const auto &f : r) {
            Add(f);
        }

        AutoReduction();

        return !r.empty();
    }

    Container data_;
};

}  // namespace groebner_basis
