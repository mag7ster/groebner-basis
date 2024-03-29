#include <cassert>
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
    explicit Monom(const std::initializer_list<Degree>& degrees_list)
        : degrees_(CreateDegrees(degrees_list)) {
    }

    auto begin() const {  // NOLINT
        return degrees_.begin();
    }

    auto end() const {  // NOLINT
        return degrees_.end();
    }

    Degree operator[](size_t index) const {
        if (index >= degrees_.size()) {
            return 0;
        }
        return degrees_[index];
    }

    bool operator==(const Monom& other) const {
        return degrees_ == other.degrees_;
    }

    bool operator!=(const Monom& other) const {
        return degrees_ != other.degrees_;
    }

    bool IsDivisibleBy(const Monom& divisor) const {

        if (CountDegrees() < divisor.CountDegrees()) {
            return false;
        }

        for (size_t i = 0; i < divisor.CountDegrees(); ++i) {
            if (divisor[i] > (*this)[i]) {
                return false;
            }
        }

        return true;
    }

    size_t CountDegrees() const {
        return degrees_.size();
    }

    Monom operator/(const Monom& other) const {

        assert(IsDivisibleBy(other));

        size_t size = CountDegrees();
        if (CountDegrees() == other.CountDegrees()) {
            size = CalcQuotientSize(other);
        }

        std::vector<Degree> newdegrees(size);
        for (size_t i = 0; i < size; ++i) {
            newdegrees[i] = (*this)[i] - other[i];
        }

        return Monom(newdegrees);
    }

    Monom operator*(const Monom& other) const {

        size_t size = std::max(CountDegrees(), other.CountDegrees());

        std::vector<Degree> newdegrees(size);
        for (size_t i = 0; i < size; ++i) {
            newdegrees[i] = (*this)[i] + other[i];
        }

        return Monom(newdegrees);
    }

private:
    size_t CalcQuotientSize(const Monom& other) const {
        size_t size = 0;
        for (size_t i = 0; i < CountDegrees(); ++i) {
            if ((*this)[i] != other[i]) {
                size = i + 1;
            }
        }
        return size;
    }

    explicit Monom(const std::vector<Degree>& vec) : degrees_(vec) {
    }

    std::vector<Degree> degrees_;
};

template <typename Stream>
Stream& operator<<(Stream& stream, const Monom& monom) {
    for (auto it = monom.begin(); it != monom.end(); ++it) {
        stream << "x" << it - monom.begin() << "^{" << *it << "}";
    }
    return stream;
}

}  // namespace groebner_basis
