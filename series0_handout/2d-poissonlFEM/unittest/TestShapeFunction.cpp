#include <gtest/gtest.h>
#include "shape.hpp"

TEST(TestShapeFunction, CornerTest) {
    // In this test, we check that the shape function
    // is 1 in each of the corners. Ie. we check that
    //
    //   lambda(i, x_j) = dirac(i, j)
    //
    // where x_j is the j-th corner of the unit triangle

    ASSERT_EQ(1.0, lambda(0, 0,
            0)) << "lambda(0, x): Should evaluate to 1 in bottom corner!";
    ASSERT_EQ(0.0, lambda(0, 1,
            0)) << "lambda(0, x): Should evaluate to 0 in right corner!";
    ASSERT_EQ(0.0, lambda(0, 0,
            1)) << "lambda(0, x): Should evaluate to 0 in top corner!";

    ASSERT_EQ(0.0, lambda(1, 0,
            0)) << "lambda(1, x): Should evaluate to 0 in bottom corner!";
    ASSERT_EQ(1.0, lambda(1, 1,
            0)) << "lambda(1, x): Should evaluate to 1 in right corner!";
    ASSERT_EQ(0.0, lambda(1, 0,
            1)) << "lambda(1, x): Should evaluate to 0 in top corner!";

    ASSERT_EQ(0.0, lambda(2, 0,
            0)) << "lambda(2, x): Should evaluate to 0 in bottom corner!";
    ASSERT_EQ(0.0, lambda(2, 1,
            0)) << "lambda(2, x): Should evaluate to 0 in right corner!";
    ASSERT_EQ(1.0, lambda(2, 0,
            1)) << "lambda(2, x): Should evaluate to 1 in top corner!";
}


TEST(TestShapeFunction, CenterTest) {
    // In this test, we check that the shape function
    // is correct in the "center" of the triangle

    ASSERT_EQ(0.25, lambda(0, 0.5,
            0.25)) << "lambda(0, 0.5, 0.25): Should evaluate to 0.25";

    ASSERT_EQ(0.5, lambda(1, 0.5,
            0.25)) << "lambda(1, 0.5, 0.25): Should evaluate to 0.5";

    ASSERT_EQ(0.25, lambda(2, 0.5,
            0.25)) << "lambda(2, 0.5, 0.25): Should evaluate to 0.25";
}

