#include <gtest/gtest.h>
#include "dirichlet_boundary.hpp"

TEST(TestDirichletBoundary, SingleTriangle) {
    // Here we check that if we only have one triangle,
    // then the every corner will be set as boundary.

    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0;

    Eigen::MatrixXi triangles(1, 3);
    triangles << 0, 1, 2;

    Eigen::VectorXi interiorIndices;

    Eigen::VectorXd u;

    u.resize(3);
    setDirichletBoundary(u, interiorIndices, vertices, triangles, [](double x, double y) {
        return x + 2 * y;
    });

    ASSERT_EQ(0, interiorIndices.size()) << "We got interior vertices when it should have been empty!";

    ASSERT_EQ(0.0, u(0)) << "did not evalute function correct in lower left corner.";
    ASSERT_EQ(1.0, u(1)) << "did not evalute function correct in right corner.";
    ASSERT_EQ(2.0, u(2)) << "did not evalute function correct in top corner.";
}

TEST(TestDirichletBoundary, FourTrianglesOneInteriorPoint) {
    // Here we check a simple configuration where we have one 
    // interior point (see "illustration" below).

    Eigen::MatrixXd vertices(5, 3);
    vertices << 0, 0, 0,
        1, -1, 0,
        1, 0, 0,
        2, 0, 0,
        1, 1, 0;


    //
    // Should look like this:
    //       4
    //      /|\
    //     / | \
    //  0 +--2--+ 3
    //     \ | /
    //      \|/
    //       1
    Eigen::MatrixXi triangles(4, 3);
    triangles << 0, 1, 2,
        1, 3, 2,
        2, 3, 4,
        0, 2, 4;

    Eigen::VectorXi interiorIndices;

    Eigen::VectorXd u;

    u.resize(5);
    setDirichletBoundary(u, interiorIndices, vertices, triangles, [](double x, double y) {
        return x + 2 * y;
    });

    ASSERT_EQ(1, interiorIndices.size()) << "We did not get the expected number of interior points!";

    ASSERT_EQ(2, interiorIndices(0)) << "Only the vertex with index 2 should be in the interior indices";

    ASSERT_EQ(0.0, u(0)) << "did not evalute function correct in lower left corner.";
    ASSERT_EQ(-1.0, u(1)) << "did not evalute function correct in right corner.";
    ASSERT_EQ(2.0, u(3)) << "did not evalute function correct in top corner.";
    ASSERT_EQ(3.0, u(4)) << "did not evalute function correct in top corner.";
}
