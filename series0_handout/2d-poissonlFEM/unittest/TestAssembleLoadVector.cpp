#include <gtest/gtest.h>
#include "load_vector.hpp"
#include "load_vector_assembly.hpp"

TEST(TestAssembleLoadVector, SingleTriangle) {
    // We test with a single triangle and make sure
    // the result is essentially the same as calling computeLoadVector

    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0;

    Eigen::MatrixXi triangles(1, 3);
    triangles << 0, 1, 2;

    Eigen::VectorXd loadVectorFromAssembly;


    // We use a slightly different f than in the square example, just to make
    // sure there is no hard-coding of functions going on!
    auto f = [](double x, double y) {
        return cos(2 * M_PI *x) * cos(2 * M_PI * y);
    };

    assembleLoadVector(loadVectorFromAssembly, vertices, triangles, f);

    Eigen::Vector3d loadVector;
    computeLoadVector(loadVector, vertices.row(0), vertices.row(1),
                      vertices.row(2), f);

    for (int i = 0; i < loadVector.rows(); ++i) {
        ASSERT_EQ(loadVector(i), loadVectorFromAssembly(i))
                << "loadVector is not the same as from assembly at component "
                << i;
    }
}

TEST(TestAssembleLoadVector, TwoTriangles) {
    // We test with two triangles, and check that the sum is computed correctly
    // the result is essentially the same as calling computeLoadVector

    Eigen::MatrixXd vertices(4, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        1, 1, 0;


    Eigen::MatrixXi triangles(2, 3);
    triangles << 0, 1, 2,
                 1, 3, 2;

    Eigen::VectorXd loadVectorFromAssembly;


    // We use a slightly different f than in the square example, just to make
    // sure there is no hard-coding of functions going on!
    auto f = [](double x, double y) {
        return cos(2 * M_PI *x) * cos(2 * M_PI * y);
    };

    assembleLoadVector(loadVectorFromAssembly, vertices, triangles, f);



    Eigen::Vector3d loadVectorLowerTriangle;
    computeLoadVector(loadVectorLowerTriangle, vertices.row(0), vertices.row(1),
                      vertices.row(2), f);

    Eigen::Vector3d loadVectorUpperTriangle;
    computeLoadVector(loadVectorUpperTriangle, vertices.row(1), vertices.row(3),
                      vertices.row(2), f);


    // In the lower left corner, only the first load vector will give the value:
    ASSERT_EQ(loadVectorLowerTriangle(0), loadVectorFromAssembly(0))
            << "Load vector not properly computed in lower left corner!";

    // In the upper right corner, only the second load vector will give the value:
    ASSERT_EQ(loadVectorUpperTriangle(1), loadVectorFromAssembly(3))
            << "Load vector not properly computed in upper right corner!";;

    // In the two values in the middle, it should be the sum of the two load
    // vectors
    ASSERT_EQ(loadVectorLowerTriangle(1) + loadVectorUpperTriangle(0),
              loadVectorFromAssembly(1))
            << "Load vector assembly does not agree on lower right corner!";

    ASSERT_EQ(loadVectorLowerTriangle(2) + loadVectorUpperTriangle(2),
              loadVectorFromAssembly(2))
            << "Load vector assembly does not agree on upper left corner!";;
}

