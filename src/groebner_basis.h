#include <cstddef>
#include <vector>
#include "functions.h"

namespace groebner_basis {

template <typename Field, typename Order = LexOrder>
class PolynomialsSet {
    using Container = std::vector<Polynom<Field, Order>>;
    using Iterator = typename Container::iterator;

public:
    PolynomialsSet(const std::initializer_list<Polynom<Field, Order>> &poly_list)
        : data_(poly_list) {
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

    void Add(const Polynom<Field, Order> &poly) {
        data_.emplace_back(poly);
    }

    void AddAt(Iterator it, const Polynom<Field, Order> &poly) {
        Add(poly);
        std::swap(*it, data_.back());
    }

    void Erase(Iterator it) {
        if (it != data_.end()) {
            std::swap(*it, data_.back());
            data_.pop_back();
        }
    }

    auto Erase(const Polynom<Field, Order> &poly) {
        auto it = std::find(data_.begin(), data_.end(), poly);
        Erase(it);
        return it;
    }

    std::optional<Polynom<Field, Order>> Reduce(const Polynom<Field, Order> &f) const {

        std::optional<Polynom<Field, Order>> temp = TryReductionForOnePass(f);
        if (!temp) {
            return std::nullopt;
        }

        Polynom<Field, Order> res{};
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
        }

        for (auto &f : (*this)) {
            f = f * Term<Field>(f.GetLargestTerm().GetCoefficient());
        }

        return count_passes;
    }

private:
    std::optional<Polynom<Field, Order>> TryReductionForOnePass(
        const Polynom<Field, Order> &f) const {

        Polynom<Field, Order> res = f;
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

    bool TryAutoReductionForOnePass() {

        auto &set = (*this);
        bool success = false;

        for (auto it = set->begin(); it != set->end(); ++it) {
            auto f = *it;
            set->Erase(it);
            auto temp = ReductionByPolynomialsSet(f, set);
            if (!temp) {
                set->AddAt(it, f);
            } else {
                success = true;
                if (temp.value() != Term<Field>(0)) {
                    set->AddAt(it, temp.value());
                }
            }
        }

        return success;
    }

    bool TryBuildGreobnerBasisForOnePass() {
        std::vector<Polynom<Field, Order>> r;
        r.reserve(Size() * (Size() - 1) / 2);
        for (auto it1 = begin(); it1 != end(); ++it1) {
            for (auto it2 = it1 + 1; it2 != end(); ++it2) {
                auto s = SPolynom(*it1, *it2);
                auto r_ij = Reduce(s);

                if (!r_ij && r_ij.value().IsZero()) {
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
