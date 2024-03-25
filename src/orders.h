#include <algorithm>
#include "monom.h"

namespace groebner_basis {
class LexOrder {
public:
    bool operator()(const Monom &a, const Monom &b) const {
        return std::lexicographical_compare(a.GetDegrees().begin(), a.GetDegrees().end(),
                                            b.GetDegrees().begin(), b.GetDegrees().end());
    }
};

}  // namespace groebner_basis
