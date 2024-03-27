#include "orders.h"

namespace groebner_basis {

template <typename Field>
class Term : public Monom {
public:
    Term(const Field& coefficient, const std::initializer_list<Degree>& degrees_list)
        : Monom(degrees_list), coef_(coefficient) {
    }

    Term(const Field& coefficient, const Monom& monom) : Monom(monom), coef_(coefficient) {
    }

    const Field& GetCoefficient() const {
        return coef_;
    }

    Term<Field> operator/(const Term<Field>& other) {
        return Term<Field>(coef_ / other.coef_,
                           static_cast<const Monom&>(*this) / static_cast<const Monom&>(other));
    }

    Term<Field> operator*(const Term<Field>& other) {
        return Term<Field>(coef_ * other.coef_,
                           static_cast<const Monom&>(*this) * static_cast<const Monom&>(other));
    }

private:
    Field coef_;
};

template <typename Field>
bool IsDeegreesEqual(const Term<Field>& m1, const Term<Field>& m2) {
    return static_cast<const Monom>(m1) == static_cast<const Monom>(m2);
}

template <typename Stream, typename Field>
Stream& operator<<(Stream& stream, const Term<Field>& term) {
    stream << term.GetCoefficient() << "*" << static_cast<const Monom&>(term);
    return stream;
}

}  // namespace groebner_basis
