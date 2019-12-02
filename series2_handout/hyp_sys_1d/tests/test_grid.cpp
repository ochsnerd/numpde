#include <gtest/gtest.h>

#include <ancse/grid.hpp>

TEST(Grid, cell_center) {
    double a = 2.0, b = 4.0;
    double dx = 2.0 / 10.0;

    int n_cells = 14;
    int n_ghost = 2;

    auto grid = Grid({a, b}, n_cells, n_ghost);

    // The cell-center of the first interior cell is
    //    a + 0.5 *dx
    ASSERT_DOUBLE_EQ(cell_center(grid, n_ghost), a + 0.5*dx);
}

TEST(Grid, PointConversions) {
  double a = -2.0, b = 4.0;

  int n_cells = 14;
  int n_ghost = 2;

  double dx = (b - a) / (n_cells - 2 * n_ghost);

  auto grid = Grid({a, b}, n_cells, n_ghost);

  double c = -1;
  double xi = reference_point(grid, c);
  // these really shouldn't be free functions, but methods of grid
  ASSERT_NEAR(c, cell_point(grid, cell_index(grid, c), xi), 1e-10);

  c = 0;
  xi = reference_point(grid, c);
  ASSERT_NEAR(c, cell_point(grid, cell_index(grid, c), xi), 1e-10);

  c = 3;
  xi = reference_point(grid, c);
  ASSERT_NEAR(c, cell_point(grid, cell_index(grid, c), xi), 1e-10);
}
