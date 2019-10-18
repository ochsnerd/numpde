#include <gtest/gtest.h>

#include <ancse/reconstruction.hpp>


// Test by checking an easy example.
TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    double ua = 1.0, ub = 2.0;
    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}

