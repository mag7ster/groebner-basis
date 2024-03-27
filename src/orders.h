#include <algorithm>
#include "monom.h"

namespace groebner_basis {
class LexOrder {
public:
    bool operator()(const Monom &a, const Monom &b) const {
        return std::lexicographical_compare(b.begin(), b.end(), a.begin(), a.end());
    }
};

}  // namespace groebner_basis
