#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ancse/reconstruction.hpp>


TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    Eigen::VectorXd ua(3), ub(3);
    ua << 1.0, 0, 1.50;
    ub << 0.1, 0, 0.15;


    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}


TEST(TestAffineReconstruction1, Example) {
  // Check with slope == 0
  auto rc = PWAffineReconstruction(make_dummy_grid(),
                                   [](Eigen::VectorXd _, double __, int ___) {return 0;});
  Eigen::MatrixXd u(2, 2);
  u << 1, 2, 3, 4;

  auto [uL, uR] = rc(u, 0);
  ASSERT_EQ(1, uL(0));
  ASSERT_EQ(2, uR(0));
  ASSERT_EQ(3, uL(1));
  ASSERT_EQ(4, uR(1));
}

TEST(TestAffineReconstruction2, Example) {
  // Check with forward slope
  auto forward_slope = [](Eigen::VectorXd u, double dx, int i)
                       {
                         return (u[i+1] - u[i])/dx;
                       };

  auto rc = PWAffineReconstruction(make_dummy_grid(), forward_slope);

  Eigen::MatrixXd u(1,3);
  u << 1,0,1;

  auto [uL,uR] = rc(u, 0);
  ASSERT_DOUBLE_EQ(uL(0), .5);
  ASSERT_DOUBLE_EQ(uR(0), -.5);
}

TEST(Testminmod_limiter, Example) {
  Eigen::VectorXd u(5);
  u << 2, 2, 1, 0, 0;

  ASSERT_DOUBLE_EQ(minmod_limiter(u, 1, 1), 0);
  ASSERT_DOUBLE_EQ(minmod_limiter(u, 1, 2), -1);
  ASSERT_DOUBLE_EQ(minmod_limiter(u, 1, 3), 0);
}
