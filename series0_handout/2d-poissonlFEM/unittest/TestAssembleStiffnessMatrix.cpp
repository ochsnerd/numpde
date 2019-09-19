#include <gtest/gtest.h>
#include "stiffness_matrix_assembly.hpp"


TEST(TestAssembleStiffnessMatrix, SingleTriangle) {
    // We test with a single triangle and make sure
    // the result is essentially the same as calling computeStiffnessMatrix

    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0;

    Eigen::MatrixXi triangles(1, 3);
    triangles << 0, 1, 2;

    SparseMatrix stiffnessMatrixFromAssembly;


    assembleStiffnessMatrix(stiffnessMatrixFromAssembly, vertices, triangles);

    Eigen::Matrix3d stiffnessMatrix;
    computeStiffnessMatrix(stiffnessMatrix, vertices.row(0), vertices.row(1),
                      vertices.row(2));

    for (int i = 0; i < stiffnessMatrix.rows(); ++i) {
        for (int j = 0; j < stiffnessMatrix.cols(); ++j) {
            ASSERT_EQ(stiffnessMatrix(i, j), stiffnessMatrixFromAssembly.coeff(i, j))
                    << "stiffnessMatrix is not the same as from assembly at component "
                    << "(" << i << ", " << j << ")";
        }
    }
}

TEST(TestAssembleStiffnessMatrix, TwoTriangles) {
    // We test with two triangles, and check that the sum is computed correctly
    // the result is essentially the same as calling computeStiffnessVector

    Eigen::MatrixXd vertices(4, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        1, 1, 0;


    Eigen::MatrixXi triangles(2, 3);
    triangles << 0, 1, 2,
                 1, 3, 2;

    SparseMatrix stiffnessMatrixFromAssembly;


    assembleStiffnessMatrix(stiffnessMatrixFromAssembly, vertices, triangles);



    Eigen::Matrix3d stiffnessMatrixLowerTriangle;
    computeStiffnessMatrix(stiffnessMatrixLowerTriangle, vertices.row(0), vertices.row(1),
                      vertices.row(2));

    Eigen::Matrix3d stiffnessMatrixUpperTriangle;
    computeStiffnessMatrix(stiffnessMatrixUpperTriangle, vertices.row(1), vertices.row(3),
                      vertices.row(2));


    // In the lower left corner, only the first load vector will give the value:
    for (int i=0; i < 3; ++i) {
        ASSERT_EQ(stiffnessMatrixLowerTriangle(0, i), stiffnessMatrixFromAssembly.coeff(0, i))
                << "stiffness matrix not properly computed in lower left corner, in combination"
                <<  "(" << 0 << ", " << i << ")";
    }

    for (int i=0; i < 3; ++i) {
        ASSERT_EQ(stiffnessMatrixLowerTriangle(i, 0), stiffnessMatrixFromAssembly.coeff(i, 0))
                << "stiffness matrix not properly computed in lower left corner, in combination"
                <<  "(" << i << ", " << 0 << ")";
    }


    // Upper corner is a bit more tricky since we need to take the ordering into
    // account
    ASSERT_EQ(stiffnessMatrixUpperTriangle(1, 0), stiffnessMatrixFromAssembly.coeff(3, 1))
                << "stiffness matrix not properly computed in upper right corner, in combination"
                <<  "(" << 3 << ", " << 1 << ")";

    ASSERT_EQ(stiffnessMatrixUpperTriangle(1, 1), stiffnessMatrixFromAssembly.coeff(3, 3))
                << "stiffness matrix not properly computed in upper right corner, in combination"
                <<  "(" << 3 << ", " << 3 << ")";

    ASSERT_EQ(stiffnessMatrixUpperTriangle(1, 2), stiffnessMatrixFromAssembly.coeff(3, 2))
                << "stiffness matrix not properly computed in upper right corner, in combination"
                <<  "(" << 3 << ", " << 2 << ")";

    ASSERT_EQ(stiffnessMatrixUpperTriangle(0, 1), stiffnessMatrixFromAssembly.coeff(1, 3))
                << "stiffness matrix not properly computed in upper right corner, in combination"
                <<  "(" << 1 << ", " << 3 << ")";

    ASSERT_EQ(stiffnessMatrixUpperTriangle(1, 1), stiffnessMatrixFromAssembly.coeff(3, 3))
                << "stiffness matrix not properly computed in upper right corner, in combination"
                <<  "(" << 3 << ", " << 3 << ")";

    ASSERT_EQ(stiffnessMatrixUpperTriangle(2, 1), stiffnessMatrixFromAssembly.coeff(2, 3))
                << "stiffness matrix not properly computed in upper right corner, in combination"
                <<  "(" << 2 << ", " << 3 << ")";

    // These should be zero:
    ASSERT_EQ(0, stiffnessMatrixFromAssembly.coeff(3, 0))
            << "StiffnessMatrix is not zero at (3,0) (but should be)!";
    ASSERT_EQ(0, stiffnessMatrixFromAssembly.coeff(0, 3))
            << "StiffnessMatrix is not zero at (0,3) (but should be)!";



    // Now comes the tricky part with the indices that are shared
    ASSERT_EQ(stiffnessMatrixLowerTriangle(1, 1)
              + stiffnessMatrixUpperTriangle(0, 0), stiffnessMatrixFromAssembly.coeff(1, 1))
            << "Wrong value at (1,1)";


    ASSERT_EQ(stiffnessMatrixLowerTriangle(1, 2)
              + stiffnessMatrixUpperTriangle(0, 2), stiffnessMatrixFromAssembly.coeff(1, 2))
            << "Wrong value at (1,2)";;


    ASSERT_EQ(stiffnessMatrixLowerTriangle(2, 1)
              + stiffnessMatrixUpperTriangle(2, 0), stiffnessMatrixFromAssembly.coeff(2, 1));

    ASSERT_EQ(stiffnessMatrixLowerTriangle(2, 2)
              + stiffnessMatrixUpperTriangle(2, 2), stiffnessMatrixFromAssembly.coeff(2, 2));

}

