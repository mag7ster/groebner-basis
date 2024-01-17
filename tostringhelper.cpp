#include "tostringhelper.h"
#include <boost/rational.hpp>

using Fraction = boost::rational<int64_t>;

namespace tostringhelper {

template <>
std::string ToString<Fraction>(Fraction x) {
    return std::to_string(x.numerator()) + "/" + std::to_string(x.denominator());
}

}  // namespace tostringhelper