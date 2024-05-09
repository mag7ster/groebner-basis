#include "orders.h"

namespace groebner_basis {

template <typename Field>
class Term : public Monom {
public:
    Term(const Field& coefficient, const std::initializer_list<Degree>& degrees_list)
        : Monom(degrees_list), coef_(coefficient) {
    }

    Term(const Field& coefficient, const Monom& monom = Monom())
        : Monom(monom), coef_(coefficient) {
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

    friend bool operator==(const Term& first, const Term& second) {
        return first.coef_ == second.coef_ && first.GetMonom() == second.GetMonom();
    }

    friend bool operator!=(const Term& first, const Term& second) {
        return !(first == second);
    }

    template <typename Stream>
    friend Stream& operator<<(Stream& stream, const Term& term) {

        if (term.GetCoefficient() == Field(1)) {
            if (term.CountSignificantDegrees() != 0) {
                stream << term.GetMonom();
            } else {
                stream << term.GetCoefficient();
            }
        } else {
            stream << term.GetCoefficient();
            if (term.CountSignificantDegrees() != 0) {
                stream << "*" << term.GetMonom();
            }
        }
        return stream;
    }

private:
    Field coef_;
};

}  // namespace groebner_basis
