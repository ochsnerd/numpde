#include <gtest/gtest.h>
#include "integrate.hpp"

TEST(IntegrationTest, ConstantIntegrationTest) {
    // In this test, we check that the integration function
    // is able to perfectly integrate a constant function

    const double constant = 4.0;
    auto f = [=](double, double) {
        return constant;
    };

    // constant * <area of triangle>
    double expectedResult = constant * 0.5;

    ASSERT_EQ(expectedResult, integrate(f)) << "Integration function does not work for constant functions!";
}

TEST(IntegrationTest, LinearIntegrationTest) {
    // In this test, we check that the integration function
    // is able to integrate a linear function

    
    auto f = [](double x, double y) {
        return 3*(x + y);
    };

    // as computed by hand:
    double expectedResult =  1.0;

    ASSERT_NEAR(expectedResult, integrate(f), 1e-8) << "Integration function does not work for linear functions!";
}

TEST(IntegrationTest, QuadraticIntegrationTest) {
    // In this test, we check that the integration function
    // is able to integrate a quadratic function


    auto f = [](double x, double y) {
        return 4 * (x*x+ 2*x*y + y*y);
    };

    // as computed by hand:
    double expectedResult = 1.0;

    ASSERT_NEAR(expectedResult, integrate(f), 1e-8) << "Integration function does not work for quadratic functions!";
}


