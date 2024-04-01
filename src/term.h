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

    const Monom& GetMonom() const {
        return *this;
    }

    Term<Field> operator/(const Term<Field>& other) const {
        return Term<Field>(coef_ / other.coef_, this->GetMonom() / other.GetMonom());
    }

    Term<Field> operator*(const Term<Field>& other) const {
        return Term<Field>(coef_ * other.coef_, this->GetMonom() * other.GetMonom());
    }

    Term<Field> operator-() const {
        return Term<Field>(-coef_, this->GetMonom());
    }

private:
    Field coef_;
};

template <typename Stream, typename Field>
Stream& operator<<(Stream& stream, const Term<Field>& term) {
    stream << term.GetCoefficient();
    if (term.CountSignificantDegrees() != 0) {
        stream << "*" << term.GetMonom();
    }
    return stream;
}

}  // namespace groebner_basis
