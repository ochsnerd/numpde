#include<gtest/gtest.h>

#include <ancse/boundary_condition.hpp>


TEST(TestBoundaryCondition, Periodic) {
    int n_cells = 10;
    int n_ghost = 2;

    Eigen::VectorXd u(n_cells);
    for(int i = 0; i < n_cells; ++i) {
        u[i] = i;
    }

    // Test you boundary condition on a small example.
    // Off-by-one errors are very common here.
}

TEST(TestBoundaryCondition, Outflow) {
    int n_cells = 10;
    int n_ghost = 2;

    Eigen::VectorXd u(n_cells);
    for(int i = 0; i < n_cells; ++i) {
        u[i] = i;
    }

    // Test you boundary condition on a small example.
    // Off-by-one errors are very common here.
}
