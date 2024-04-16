#include <initializer_list>
#include <utility>
#include <vector>
#include "polynom.h"

namespace groebner_basis {

template <typename Field, typename Order = LexOrder>
class PolynomialsSet {
public:
    using Container = std::vector<Polynom<Field, Order>>;
    using Iterator = typename Container::iterator;

    PolynomialsSet(const std::initializer_list<Polynom<Field, Order>>& poly_list)
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

    void Add(const Polynom<Field, Order>& poly) {
        data_.emplace_back(poly);
    }

    void AddAt(Iterator it, const Polynom<Field, Order>& poly) {
        Add(poly);
        std::swap(*it, data_.back());
    }

    void Erase(Iterator it) {
        if (it != data_.end()) {
            std::swap(*it, data_.back());
            data_.pop_back();
        }
    }

    auto Erase(const Polynom<Field, Order>& poly) {
        auto it = std::find(data_.begin(), data_.end(), poly);
        Erase(it);
        return it;
    }

private:
    Container data_;
};

}  // namespace groebner_basis
