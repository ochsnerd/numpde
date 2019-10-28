#include <Eigen/Dense>
#include <ancse/cfl_condition.hpp>
#include <gtest/gtest.h>

TEST(TestCFL, Example) {
  auto cfl1 = ConcreteCFL(make_dummy_grid(),
                          make_dummy_model(),
                          1);
  auto cfl2 = ConcreteCFL(make_dummy_grid(),
                          make_dummy_model(),
                          2);

  ASSERT_DOUBLE_EQ(cfl1(Eigen::Vector2d{0.01, 0}), 1);
  ASSERT_DOUBLE_EQ(cfl1(Eigen::Vector2d{-0.01, 0}), 1);
  ASSERT_DOUBLE_EQ(cfl2(Eigen::Vector2d{0.01, 0}), 2);
}
