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

TEST(TestEulerModel2, Conversions) {
  auto model = Euler();
  Eigen::VectorXd u;
  u = Eigen::Vector3d(1,1,1);

  ASSERT_DOUBLE_EQ(u(0), model.cons_to_prim(model.prim_to_cons(u))(0));
  ASSERT_DOUBLE_EQ(u(1), model.cons_to_prim(model.prim_to_cons(u))(1));
  ASSERT_DOUBLE_EQ(u(2), model.cons_to_prim(model.prim_to_cons(u))(2));

  u = Eigen::Vector3d(1,0,0);

  ASSERT_DOUBLE_EQ(u(0), model.cons_to_prim(model.prim_to_cons(u))(0));
  ASSERT_DOUBLE_EQ(u(1), model.cons_to_prim(model.prim_to_cons(u))(1));
  ASSERT_DOUBLE_EQ(u(2), model.cons_to_prim(model.prim_to_cons(u))(2));

  u = Eigen::Vector3d(M_PI, -12.34, 42);

  ASSERT_DOUBLE_EQ(u(0), model.cons_to_prim(model.prim_to_cons(u))(0));
  ASSERT_DOUBLE_EQ(u(1), model.cons_to_prim(model.prim_to_cons(u))(1));
  ASSERT_DOUBLE_EQ(u(2), model.cons_to_prim(model.prim_to_cons(u))(2));

  ASSERT_DOUBLE_EQ(model.cons_to_prim(u)(0), u(0));
  ASSERT_DOUBLE_EQ(model.cons_to_prim(u)(1), u(1) / u(0));
  ASSERT_DOUBLE_EQ(model.cons_to_prim(u)(2), (u(2) - .5 * u(1) * u(1) / u(0)) * (model.get_gamma() - 1));

  ASSERT_DOUBLE_EQ(model.prim_to_cons(u)(0), u(0));
  ASSERT_DOUBLE_EQ(model.prim_to_cons(u)(1), u(0) * u(1));
  ASSERT_DOUBLE_EQ(model.prim_to_cons(u)(2), u(2) / (model.get_gamma() - 1) + 0.5 * u(0) * u(1) * u(1));
}
