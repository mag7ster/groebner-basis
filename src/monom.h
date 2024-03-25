#include <cstddef>
#include <vector>
#include <initializer_list>

namespace groebner_basis {

using Degree = unsigned int;

static std::vector<Degree> CreateDegrees(const std::initializer_list<Degree>& degrees_list) {

    size_t size = 0, zeros = 0;

    for (const auto& elem : degrees_list) {
        if (elem == 0) {
            ++zeros;
        } else {
            size += zeros + 1;
            zeros = 0;
        }
    }

    std::vector<Degree> degrees(size);
    for (auto it = degrees_list.begin(); it != degrees_list.end(); ++it) {
        degrees[it - degrees_list.begin()] = *it;
    }

    return degrees;
}

class Monom {
public:
    Monom(const std::initializer_list<Degree>& degrees_list)
        : degrees_(CreateDegrees(degrees_list)) {
    }

    const std::vector<Degree>& GetDegrees() const {
        return degrees_;
    }

    bool IsDivisibleBy(const Monom& divisor) const {
        for (size_t i = 0; i < std::min(degrees_.size(), divisor.degrees_.size()); ++i) {
            if (divisor.degrees_[i] > degrees_[i]) {
                return false;
            }
        }
        return true;
    }

private:
    std::vector<Degree> degrees_;
};

template <typename Stream>
Stream& operator<<(Stream& stream, const Monom& monom) {
    for (auto it = monom.GetDegrees().begin(); it != monom.GetDegrees().end(); ++it) {
        stream << "x" << it - monom.GetDegrees().begin() << "^{" << *it << "}";
    }
    return stream;
}

bool IsDeegreesEqual(const Monom& m1, const Monom& m2) {
    return m1.GetDegrees() == m2.GetDegrees();
}

}  // namespace groebner_basis
