#include <gtest/gtest.h>

#include "load_vector.hpp"

TEST(TestElementVector, TestConstantF) {
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0;

    const double constant = 1.0;
    auto f = [=](double, double) {
        return constant;
    };

    Eigen::Vector3d elementVector;

    computeLoadVector(elementVector, vertices.row(0), vertices.row(1), vertices.row(2), f);

    // This should be equal to
    // 
    //   int_K lambda(0,x) * f(x) dx
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected0 = constant / 6.0;
    ASSERT_NEAR(expected0, elementVector(0),  1e-8) << "Element vector not correct for lambda(0)";

    // This should be equal to
    // 
    //   int_K lambda(1,x) * f(x) dx 
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected1 = constant / 6.0;
    ASSERT_NEAR(expected1, elementVector(1), 1e-8) << "Element vector not correct for lambda(1)";

    // This should be equal to
    // 
    //   int_K lambda(2,x) * f(x) dx
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected2 = constant / 6.0;
    ASSERT_NEAR(expected2, elementVector(2), 1e-8) << "Element vector not correct for lambda(2)";
}


TEST(TestElementVector, TestConstantFSkewTriangle) {
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0.4, 0,
        0, 1, 0;

    const double constant = 1.0;
    auto f = [=](double, double) {
        return constant;
    };

    Eigen::Vector3d elementVector;

    computeLoadVector(elementVector, vertices.row(0), vertices.row(1), vertices.row(2), f);

    // This should be equal to
    // 
    //   int_K lambda(0,x) * f(x) dx = 0.5 * 1.0 /6.0
    // 
    // where K is the triangle spanned by (0,0), (1, 0.4), (1,0)
    const double expected0 =  constant / 6;
    ASSERT_NEAR(expected0, elementVector(0), 1e-8);

    // This should be equal to
    // 
    //   int_K lambda(1,x) * f(x) dx 
    // 
    // where K is the triangle spanned by (0,0), (1, 0.4), (1,0)
    const double expected1 =  constant / 6.;
    ASSERT_NEAR(expected1, elementVector(1), 1e-8);

    // This should be equal to
    // 
    //   int_K lambda(2,x) * f(x) dx 
    // 
    // where K is the triangle spanned by (0,0), (1, 0.4), (1,0)
    const double expected2 = constant / 6.;
    ASSERT_NEAR(expected2, elementVector(2), 1e-8);
}

TEST(TestElementVector, TestLinearFSkewTriangle) {
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0.4, 0,
        0, 1, 0;

    
    auto f = [=](double x, double y) {
        return x + y;
    };

    Eigen::Vector3d elementVector;

    computeLoadVector(elementVector, vertices.row(0), vertices.row(1), vertices.row(2), f);

    // This should be equal to
    // 
    //   int_K lambda(0,x) * f(x) dx = 0.5 * 1.0 /6.0
    // 
    // where K is the triangle spanned by (0,0), (1, 0.4), (1,0)
    const double expected0 = 0.1;
    ASSERT_NEAR(expected0, elementVector(0), 1e-8);

    // This should be equal to
    // 
    //   int_K lambda(1,x) * f(x) dx
    // 
    // where K is the triangle spanned by (0,0), (1, 0.4), (1,0)
    const double expected1 = 0.158333333333333333333;
    ASSERT_NEAR(expected1, elementVector(1), 1e-8);

    // This should be equal to
    // 
    //   int_K lambda(2,x) * f(x) dx 
    // 
    // where K is the triangle spanned by (0,0), (1, 0.4), (1,0)
    const double expected2 = 0.14166666666666667;
    ASSERT_NEAR(expected2, elementVector(2), 1e-8);
}