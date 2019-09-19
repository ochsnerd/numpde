#include <gtest/gtest.h>
#include "coordinate_transform.hpp"

TEST(TestCoordinateTransform, IdentityTest) {
    // We should make sure the coordinate transform 
    // gives the identity when we supply it with the
    // standard basis {(1,0), (0, 1)}.

    Eigen::Vector3d e1(1, 0, 0);
    Eigen::Vector3d e2(0, 1, 0);

    auto transform = makeCoordinateTransform(e1, e2);

    ASSERT_EQ(1, transform(0, 0));
    ASSERT_EQ(0, transform(0, 1));
    ASSERT_EQ(0, transform(1, 0));
    ASSERT_EQ(1, transform(1, 1));
}

TEST(TestCoordinateTransform, SkewTest) {
    // We should make sure the coordinate transform 
    // gives the correct version when we supply it with 
    // non-straight triangle

    Eigen::Vector3d a(0, 0, 0);
    Eigen::Vector3d b(1, 0.4, 0);
    Eigen::Vector3d c(0, 1, 0);

    auto transform = makeCoordinateTransform(b - a, c - a);

    ASSERT_EQ(1, transform(0, 0));
    ASSERT_EQ(0, transform(0, 1));
    ASSERT_EQ(0.4, transform(1, 0));
    ASSERT_EQ(1, transform(1, 1));

    Eigen::Vector2d e1Transformed = transform * Eigen::Vector2d(1.0, 0.0);

    ASSERT_EQ(1, e1Transformed(0)) ;
    ASSERT_NEAR(0.4, e1Transformed(1), 1e-8);

    Eigen::Vector2d e2Transformed = transform * Eigen::Vector2d(0.0, 1.0);

    ASSERT_EQ(0, e2Transformed(0));
    ASSERT_EQ(1, e2Transformed(1));

}
