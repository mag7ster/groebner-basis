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

private:
    Field coef_;
};

template <typename Stream, typename Field>
Stream& operator<<(Stream& stream, const Term<Field>& term) {
    stream << term.GetCoefficient() << "*" << static_cast<Monom>(term);
    return stream;
}

}  // namespace groebner_basis
