#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <ancse/model.hpp>

TEST(TestEulerModel, Example) {
  auto model = Euler();
  auto u = Eigen::Vector3d(1,1,1);
  auto gamma = model.get_gamma();
  auto c = std::sqrt(0.5*gamma*(gamma-1));
  auto p = 0.5 * (gamma - 1);
  auto H = 0.5 * (gamma + 1);

  // Make sure we get the right values for u = (1,1,1)
  //Fluxes
  ASSERT_DOUBLE_EQ(model.flux(u)(0), 1.);
  ASSERT_DOUBLE_EQ(model.flux(u)(1), 1. + p);
  ASSERT_DOUBLE_EQ(model.flux(u)(2), 1. + p);

  // Eigenvalues
  ASSERT_DOUBLE_EQ(model.eigenvalues(u)(0), 1. - c);
  ASSERT_DOUBLE_EQ(model.eigenvalues(u)(1), 1.);
  ASSERT_DOUBLE_EQ(model.eigenvalues(u)(2), 1. + c);

  // Eigenvectors
  auto evs = Eigen::Matrix3d();
  evs <<
    1.,     1.,   1.,
    1. - c, 1.,   1. + c,
    H - c,  0.5,  H + c;
  ASSERT_DOUBLE_EQ((model.eigenvectors(u) - evs).norm(), 0);
}
