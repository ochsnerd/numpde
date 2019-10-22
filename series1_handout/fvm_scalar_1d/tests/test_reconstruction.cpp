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

TEST(Testminmod_limiter, Example) {
  Eigen::VectorXd u(5);
  u << 2, 2, 1, 0, 0;

  ASSERT_DOUBLE_EQ(minmod_limiter(u, 1, 1), 0);
  ASSERT_DOUBLE_EQ(minmod_limiter(u, 1, 2), -1);
  ASSERT_DOUBLE_EQ(minmod_limiter(u, 1, 3), 0);
}

TEST(Testsuperbee_limiter, Example) {
  ASSERT_EQ("Did you find some test-examples for superbee?","no :(");
}

// Test mathematical functions
TEST(Testslope, Example) {
  ASSERT_DOUBLE_EQ(slope(1,1,1), 0);
  ASSERT_DOUBLE_EQ(slope(0,1,1), 1);
  ASSERT_DOUBLE_EQ(slope(1,0,1), -1);
  ASSERT_DOUBLE_EQ(slope(-2,2,4), 1);
  ASSERT_DOUBLE_EQ(slope(0.0,M_PI,M_PI), 1);
}

TEST(Testminmod, Example) {
  ASSERT_EQ(minmod(-1,1), 0);
  ASSERT_EQ(minmod(1,1), 1);
  ASSERT_EQ(minmod(-1,-1), -1);
  ASSERT_EQ(minmod(-1,-2), -1);
  // depends on definiton of sign(0)
  // or does it? maybe just for maxmod
  ASSERT_EQ(minmod(1,0), 0);
  ASSERT_EQ(minmod(-1, 0), 0);
}

TEST(Testmaxmod, Example) {
  ASSERT_EQ(maxmod(-1,1), 0);
  ASSERT_EQ(maxmod(1,1), 1);
  ASSERT_EQ(maxmod(-1,-1), -1);
  ASSERT_EQ(maxmod(-1,-2), -2);
  // depends on definiton of sign(0)
  ASSERT_EQ(maxmod(1,0), 1);
  ASSERT_EQ(maxmod(-1, 0), 0);
}
