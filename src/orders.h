#include <algorithm>
#include <numeric>
#include "monom.h"

namespace groebner_basis {
class LexOrder {
public:
    bool operator()(const Monom &a, const Monom &b) const {
        return std::lexicographical_compare(b.begin(), b.end(), a.begin(), a.end());
    }
};

class GrevLexOrder {
public:
    bool operator()(const Monom &a, const Monom &b) const {
        Monom::Degree sum1 = std::accumulate(a.begin(), a.end(), 0),
                      sum2 = std::accumulate(b.begin(), b.end(), 0);
        if (sum1 == sum2) {
            return LexOrder()(a, b);
        }
        return sum1 > sum2;
    }
};

}  // namespace groebner_basis
