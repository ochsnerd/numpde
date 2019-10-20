#include <gtest/gtest.h>

#include <ancse/reconstruction.hpp>


// Test by checking an easy example.
TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    double ua = 1.0, ub = 2.0;
    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}

TEST(TestSecondOrderReconstruction, Example) {
  // Check with slope == 0
  auto rc = SecondOrderReconstruction(make_dummy_grid(),
                                       [](Eigen::VectorXd _, double __, int ___) {return 0;});
  Eigen::VectorXd u(2);
  u << 1, 2;

  auto [uL, uR] = rc(u, 0);
  ASSERT_EQ(1, uL);
  ASSERT_EQ(2, uR);
}

TEST(SecondOrderReconstruction, Example) {
  // Check with forward slope
  auto forward_slope = [](Eigen::VectorXd u, double dx, int i)
                       {
                         return (u[i+1] - u[i])/dx;
                       };
  auto rc = SecondOrderReconstruction(make_dummy_grid(), forward_slope);

  Eigen::VectorXd u(3);
  u << 1,0,1;

  auto [uL,uR] = rc(u, 0);
  ASSERT_DOUBLE_EQ(uL, .5);
  ASSERT_DOUBLE_EQ(uR, -.5);
}
