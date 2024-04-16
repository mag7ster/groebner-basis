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

    Term operator/(const Term& other) const {
        return Term(coef_ / other.coef_, this->GetMonom() / other.GetMonom());
    }

    Term operator*(const Term& other) const {
        return Term(coef_ * other.coef_, this->GetMonom() * other.GetMonom());
    }

    Term operator-() const {
        return Term(-coef_, this->GetMonom());
    }

    bool operator==(const Term& other) {
        return coef_ == other.coef_ && GetMonom() == other.GetMonom();
    }

    bool operator!=(const Term& other) {
        return !(*this == other);
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
