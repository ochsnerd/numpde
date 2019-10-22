#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include <ancse/reconstruction.hpp>

template <class SlopeLimiter>
void check_TVD(const SlopeLimiter& sigma,
               const Eigen::VectorXd& u,
               std::pair<double,double> bounds,
               double dx) {
  for (int i = 1; i < u.size() - 1; ++i) {
    auto jump = dx * (sigma(u, dx, i)) / (u[i+1] - u[i]);

    EXPECT_LE(jump, bounds.second);
    EXPECT_GE(jump, bounds.first);
  }
}

template <class SlopeLimiter>
void check_TVD_sin(const SlopeLimiter& sigma, std::pair<double, double> bounds) {
  int n = 1000;
  double dx = 0.01;

  Eigen::VectorXd u(n);

  for (int i = 0; i < u.size(); ++i) {
    u[i] = std::sin(2 * M_PI * i / n);
  }

  check_TVD(sigma, u, bounds, dx);
}

template <class SlopeLimiter>
void check_TVD_sin_perturbed(const SlopeLimiter& sigma, std::pair<double, double> bounds) {
  int n = 1000;
  double dx = 0.01;

  Eigen::VectorXd u(n);

  std::mt19937 eng(42);
  std::normal_distribution<double> N(0, 0.1);

  for (int i = 0; i < u.size(); ++i) {
    u[i] = std::sin(2 * M_PI * i / n) + N(eng);
  }

  check_TVD(sigma, u, bounds, dx);
}

template <class SlopeLimiter>
void check_TVD_local_extrema(const SlopeLimiter& sigma) {
  Eigen::VectorXd maximum(3);
  maximum << 0, 1, 0;

  Eigen::VectorXd minimum(3);
  minimum << 0, -1, 0;

  EXPECT_DOUBLE_EQ(sigma(maximum, 1, 1), 0);
  EXPECT_DOUBLE_EQ(sigma(minimum, 1, 1), 0);
}

template <class SlopeLimiter>
void check_TVD(const SlopeLimiter& sigma) {
  check_TVD_sin(sigma, {-2, 2});
  check_TVD_sin_perturbed(sigma, {-2, 2});
  check_TVD_local_extrema(sigma);
}

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
