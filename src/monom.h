#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <vector>
#include <initializer_list>
#include <type_traits>

namespace groebner_basis {

template <typename Iterator, typename T>
concept IsIteratorValueEqualsT =
    std::is_same_v<T, typename std::iterator_traits<Iterator>::value_type>;

class Monom {
public:
    using Degree = unsigned int;

    Monom() = default;

    explicit Monom(std::initializer_list<Degree> degrees_list)
        : Monom(TrimTrailingZeros(degrees_list.begin(), degrees_list.end())) {
    }

    static Monom BuildFromVectorDegrees(const std::vector<Monom::Degree>& vector_degrees) {
        return Monom(Monom::TrimTrailingZeros(vector_degrees.begin(), vector_degrees.end()));
    }

    auto begin() const {  // NOLINT
        return degrees_->begin();
    }

    auto end() const {  // NOLINT
        return degrees_->end();
    }

    auto rbegin() const {  // NOLINT
        return degrees_->rbegin();
    }

    auto rend() const {  // NOLINT
        return degrees_->rend();
    }

    Degree operator[](size_t index) const {
        if (index >= degrees_->size()) {
            return 0;
        }
        return degrees_->at(index);
    }  // Я не знаю, как сделать это по другому. Если я буду возвращать const vector&, то не смогу
       // заифать, чтобы возвращался 0 когда степень больше чем есть в векторе

    bool operator==(const Monom& other) const {
        return *degrees_ == *other.degrees_;
    }

    bool operator!=(const Monom& other) const {
        return !(*this == other);
    }

    bool IsDivisibleBy(const Monom& divisor) const {

        if (FirstIndexAfterLastNonZeroDegree() < divisor.FirstIndexAfterLastNonZeroDegree()) {
            return false;
        }

        return std::mismatch(divisor.begin(), divisor.end(), begin(),
                             [](Degree divisor_degree, Degree divident_degree) {
                                 return divisor_degree <= divident_degree;
                             })
                   .first == divisor.end();
    }

    size_t FirstIndexAfterLastNonZeroDegree() const {
        return degrees_->size();
    }

    Monom operator/(const Monom& other) const {

        assert(IsDivisibleBy(other));

        size_t size = FirstIndexAfterLastNonZeroDegree();
        if (FirstIndexAfterLastNonZeroDegree() == other.FirstIndexAfterLastNonZeroDegree()) {
            size = CalcQuotientSize(other);
        }

        std::vector<Degree> newdegrees(size);
        for (size_t i = 0; i < size; ++i) {
            newdegrees[i] = (*this)[i] - other[i];
        }

        return Monom(std::move(newdegrees));
    }

    Monom operator*(const Monom& other) const {

        size_t size =
            std::max(FirstIndexAfterLastNonZeroDegree(), other.FirstIndexAfterLastNonZeroDegree());

        std::vector<Degree> newdegrees(size);
        for (size_t i = 0; i < size; ++i) {
            newdegrees[i] = (*this)[i] + other[i];
        }

        return Monom(std::move(newdegrees));
    }

    friend Monom LCM(const Monom& m1, const Monom& m2);

private:
    size_t CalcQuotientSize(const Monom& other) const {

        size_t size = std::mismatch(degrees_->rbegin(), degrees_->rend(), other.degrees_->rbegin(),
                                    other.degrees_->rend())
                          .first.base() -
                      begin();
        return size;
    }

    explicit Monom(std::vector<Degree>&& vec)
        : degrees_(std::make_shared<std::vector<Degree>>(std::move(vec))) {
    }

    template <typename Iterator>
        requires IsIteratorValueEqualsT<Iterator, Degree>
    static std::vector<Degree> TrimTrailingZeros(Iterator begin, Iterator end) {

        std::vector<Degree> zero = {0};

        auto trailingzeros =
            std::find_end(begin, end, zero.begin(), zero.end(), std::not_equal_to<Degree>());

        if (trailingzeros == end) {
            trailingzeros = begin;
        } else {
            ++trailingzeros;
        }

        std::vector<Degree> degrees(std::distance(begin, trailingzeros));

        std::copy(begin, trailingzeros, degrees.begin());
        return degrees;
    }

    std::shared_ptr<const std::vector<Degree>> degrees_ =
        std::make_shared<const std::vector<Degree>>();
};

template <typename Stream>
Stream& operator<<(Stream& stream, const Monom& monom) {

    for (auto it = monom.begin(); it != monom.end(); ++it) {
        if (*it) {
            stream << "x" << it - monom.begin();
            if (*it - 1) {
                stream << "^{" << *it << "}";
            }
        }
    }

    return stream;
}

}  // namespace groebner_basis
