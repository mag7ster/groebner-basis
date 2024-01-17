#include "groebner_basis.h"

namespace groebner_basis {

Monom::Monom(std::initializer_list<int> list) : degs_(list) {
    while (!degs_.empty() && degs_.back() == 0) {
        degs_.pop_back();
    }
}

bool LexOrder::operator()(const Monom &a, const Monom &b) {
    for (size_t i = 0; i < std::min(a.degs_.size(), b.degs_.size()); ++i) {
        if (a.degs_[i] == b.degs_[i]) {
            continue;
        }
        return a.degs_[i] > b.degs_[i];
    }
    return a.degs_.size() > b.degs_.size();
}

}  // namespace groebner_basis