#include <gtest/gtest.h>
#include "grad_shape.hpp"

TEST(TestGradientShapeFunction, ValueTest) {
    // We check that the gradient is correct for each shape function
    // We know that
    //
    //   grad Lambda_0 = (-1, -1)
    //
    //   grad Lambda_1 = (1,   0)
    //
    //   grad Lambda_2 = (0,   1)
    //
    // and it should be constant


    {
        auto gradLambda0 = gradientLambda(0, 0, 0);
        ASSERT_EQ(-1, gradLambda0(0)) << "gradLambda(0, x) is not -1 in first component";
        ASSERT_EQ(-1, gradLambda0(1)) << "gradLambda(0, x) is not -1 in second component";
    }

    {
        auto gradLambda1 = gradientLambda(1, 0, 0);
        ASSERT_EQ(1, gradLambda1(0)) << "gradLambda(1, x) is not 1 in first component";
        ASSERT_EQ(0, gradLambda1(1)) << "gradLambda(1, x) is not 0 in second component";
    }

    {
        auto gradLambda2 = gradientLambda(2, 0, 0);
        ASSERT_EQ(0, gradLambda2(0)) << "gradLambda(2, x) is not 0 in first component";
        ASSERT_EQ(1, gradLambda2(1)) << "gradLambda(2, x) is not 1 in second component";
    }
}
