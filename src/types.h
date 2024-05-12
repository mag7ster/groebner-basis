
#include <cassert>

namespace groebner_basis {

template <typename T, T Tmod>
class Modulus {
public:
    Modulus() : value_(0) {
    }

    Modulus(T value) : value_(value % Tmod) {
        if (value_ < 0) {
            value_ += Tmod;
        }
    }

    Modulus operator-() const {
        return Modulus(Tmod - value_);
    }

    Modulus& operator+=(Modulus other) {
        return (*this) = (*this) + other;
    }

    Modulus& operator-=(Modulus other) {
        return (*this) = (*this) - other;
    }

    Modulus& operator*=(Modulus other) {
        return (*this) = (*this) * other;
    }

    Modulus& operator/=(Modulus other) {
        return (*this) = (*this) / other;
    }

    // operator bool() const {
    //     return value_;
    // }

    friend Modulus operator+(Modulus first, Modulus second) {
        return Modulus((first.value_ + second.value_) % Tmod);
    }

    friend Modulus operator-(Modulus first, Modulus second) {
        return first + (-second);
    }

    friend Modulus operator*(Modulus first, Modulus second) {
        return Modulus((first.value_ * second.value_) % Tmod);
    }

    friend Modulus operator/(Modulus first, Modulus second) {
        return first * second.Inverse();
    }

    friend bool operator==(Modulus first, Modulus second) {
        return first.value_ == second.value_;
    }

    friend bool operator!=(Modulus first, Modulus second) {
        return !(first == second);
    }

    friend bool operator<(Modulus first, Modulus second) {
        return first.value_ < second.value_;
    }

    friend bool operator<=(Modulus first, Modulus second) {
        return first < second || first == second;
    }

    friend bool operator>(Modulus first, Modulus second) {
        return first >= second && first != second;
    }

    friend bool operator>=(Modulus first, Modulus second) {
        return !(first < second);
    }

    template <typename Stream>
    friend Stream& operator<<(Stream& stream, Modulus modulus) {
        stream << modulus.value_;
        return stream;
    }

private:
    Modulus BinPow(T degree) const {
        Modulus result = 1, current = (*this);
        while (degree != 0) {
            if (degree % 2 == 1) {
                result *= current;
            }
            current *= current;
            degree /= 2;
        }
        return result;
    }

    Modulus Inverse() const {
        assert(value_);
        return BinPow(Tmod - 2);
    }

    T value_;
};
}  // namespace groebner_basis
