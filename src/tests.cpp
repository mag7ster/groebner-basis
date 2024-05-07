#include "groebner_basis.h"
#include <boost/rational.hpp>
#include <gtest/gtest.h>

using Fraction = boost::rational<int64_t>;

TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}
