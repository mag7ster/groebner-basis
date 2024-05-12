#include "groebner_basis.h"
#include "types.h"

#include <gtest/gtest.h>
#include <boost/rational.hpp>
#include <fstream>

namespace {

namespace gb = groebner_basis;

using ModInt = gb::Modulus<std::int64_t, 998244353>;
using Fraction = boost::rational<std::int64_t>;

void Check(gb::PolynomialsSet<ModInt>& find, gb::PolynomialsSet<ModInt>& ans) {
    find.BuildGreobnerBasis();
    std::sort(ans.begin(), ans.end());
    EXPECT_EQ(find, ans);
}

void CheckFromFile() {
    int state;
    std::string str;

    gb::PolynomialsSet<ModInt> find, ans;
    std::ifstream file("../tests.txt");

    EXPECT_EQ(file.is_open(), true);

    int num = 0;

    while (file >> str) {

        if (str == "test") {
            state = 0;
            continue;
        }
        if (str == "Groebner") {
            state = 1;
            continue;
        }
        if (str == "calc") {
            Check(find, ans);
            find.Clear();
            ans.Clear();
            continue;
        }

        gb::Polynom<ModInt> poly = gb::Polynom<ModInt>::BuildFromString(str);

        if (state == 0) {
            find.Add(poly);
        } else {
            ans.Add(poly);
        }
    }
}
}  // namespace

TEST(GroebnerBasisTest, Stress) {
    CheckFromFile();
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}
