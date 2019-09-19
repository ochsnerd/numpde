#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "stiffness_matrix.hpp"

TEST(TestStiffnessMatrix, TestReferenceElement) {
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0;
    
    Eigen::Matrix3d stiffnessMatrix;


    computeStiffnessMatrix(stiffnessMatrix, vertices.row(0), vertices.row(1), vertices.row(2));

    // This should be equal to
    // 
    //   int_K grad lambda(0, x) * grad lambda(0, x) dx = 1
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected00 = 1;
    ASSERT_NEAR(expected00, stiffnessMatrix(0, 0), 1e-8) << "stiffness matrix vector not correct at (0,0)";


    // This should be equal to
    // 
    //   int_K grad lambda(0, x) * grad lambda(1, x) dx = -0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected01 = -0.5;
    ASSERT_NEAR(expected01, stiffnessMatrix(0, 1), 1e-8) << "stiffness matrix vector not correct at (0,1)";

    // This should be equal to
    // 
    //   int_K grad lambda(0, x) * grad lambda(2, x) dx = -0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected02 = -0.5;
    ASSERT_NEAR(expected02, stiffnessMatrix(0, 2), 1e-8) << "stiffness matrix vector not correct at (0,2)";

    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(0, x) dx = -0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected10 = -0.5;
    ASSERT_NEAR(expected10, stiffnessMatrix(1, 0), 1e-8) << "stiffness matrix vector not correct at (1,0)";


    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(1, x) dx = 0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected11 = 0.5;
    ASSERT_NEAR(expected11, stiffnessMatrix(1, 1), 1e-8) << "stiffness matrix vector not correct at (1,1)";

    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(2, x) dx = 0
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected12 = 0;
    ASSERT_NEAR(expected12, stiffnessMatrix(1, 2), 1e-8) << "stiffness matrix vector not correct at (1,2)";


    // This should be equal to
    // 
    //   int_K grad lambda(2, x) * grad lambda(0, x) dx = -0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected20 = -0.5;
    ASSERT_NEAR(expected20, stiffnessMatrix(2, 0), 1e-8) << "stiffness matrix vector not correct at (2, 0)";


    // This should be equal to
    // 
    //   int_K grad lambda(2, x) * grad lambda(1, x) dx = 0
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected21 = 0;
    ASSERT_NEAR(expected21, stiffnessMatrix(2, 1), 1e-8) << "stiffness matrix vector not correct at (2,1)";

    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(2, x) dx = 0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0), (1,0)
    const double expected22 = 0.5;
    ASSERT_NEAR(expected22, stiffnessMatrix(2, 2), 1e-8) << "stiffness matrix vector not correct at (2,2)";
}


TEST(TestStiffnessMatrix, TestSkewTriagnle) {
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0.4, 0,
        0, 1, 0;

    Eigen::Matrix3d stiffnessMatrix;

    computeStiffnessMatrix(stiffnessMatrix, vertices.row(0), vertices.row(1), vertices.row(2));

    // This should be equal to
    // 
    //   int_K grad lambda(0, x) * grad lambda(0, x) dx = 0.68
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected00 = 0.68;
    ASSERT_NEAR(expected00, stiffnessMatrix(0, 0), 1e-8) << "stiffness matrix vector not correct at (0,0)";


    // This should be equal to
    // 
    //   int_K grad lambda(0, x) * grad lambda(1, x) dx = -0.3
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected01 = -0.3;
    ASSERT_NEAR(expected01, stiffnessMatrix(0, 1), 1e-8) << "stiffness matrix vector not correct at (0,1)";

    // This should be equal to
    // 
    //   int_K grad lambda(0, x) * grad lambda(2, x) dx = -0.38
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected02 = -0.38;
    ASSERT_NEAR(expected02, stiffnessMatrix(0, 2), 1e-8) << "stiffness matrix vector not correct at (0,2)";

    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(0, x) dx = -0.3
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected10 = -0.3;
    ASSERT_NEAR(expected10, stiffnessMatrix(1, 0), 1e-8) << "stiffness matrix vector not correct at (1,0)";


    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(1, x) dx = 0.5
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected11 = 0.5;
    ASSERT_NEAR(expected11, stiffnessMatrix(1, 1), 1e-8) << "stiffness matrix vector not correct at (1,1)";

    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(2, x) dx = 0.2
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected12 = -0.2;
    ASSERT_NEAR(expected12, stiffnessMatrix(1, 2), 1e-8) << "stiffness matrix vector not correct at (1,2)";


    // This should be equal to
    // 
    //   int_K grad lambda(2, x) * grad lambda(0, x) dx = -0.38
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected20 = -0.38;
    ASSERT_NEAR(expected20, stiffnessMatrix(2, 0), 1e-8) << "stiffness matrix vector not correct at (2, 0)";


    // This should be equal to
    // 
    //   int_K grad lambda(2, x) * grad lambda(1, x) dx = -0.2
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected21 = -0.2;
    ASSERT_NEAR(expected21, stiffnessMatrix(2, 1), 1e-8) << "stiffness matrix vector not correct at (2,1)";

    // This should be equal to
    // 
    //   int_K grad lambda(1, x) * grad lambda(2, x) dx = 0.58
    // 
    // where K is the triangle spanned by (0,0), (1,0.4), (1,0)
    const double expected22 = 0.58;
    ASSERT_NEAR(expected22, stiffnessMatrix(2, 2), 1e-8) << "stiffness matrix vector not correct at (2,2)";
}
