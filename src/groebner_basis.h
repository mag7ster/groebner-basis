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

    void AutoReduction() {
        while (TryAutoReductionForOnePass()) {
            // pass
        }
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

    Container data_;
};

}  // namespace groebner_basis
