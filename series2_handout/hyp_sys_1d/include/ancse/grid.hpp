#ifndef HYPSYS1D_GRID_HPP
#define HYPSYS1D_GRID_HPP

#include <utility>
#include <cmath>

/// A uniform Cartesian Grid in 1D.
/** The physical domain is [a, b], where `auto [a, b] = domain;`. The
 *  ghost-cells, if present, lie outside the physical domain.
 */
struct Grid {
    Grid(std::pair<double, double> domain, int n_cells, int n_ghost);

    std::pair<double, double> domain;
    int n_cells; ///< number of cells, including ghost-cells.
    int n_ghost; ///< number of ghost_cells on either side.
    double dx;   ///< mesh width of a uniform grid
};

/// Compute the center of cell `i`.
double cell_center(const Grid &grid, int i);

/// Compute the physical point of cell `i` from given reference point.
double cell_point(const Grid &grid, int i, double xi);

/// Compute the reference point in [0,1] given a physical point.
double reference_point(const Grid& grid, double x);

/// Compute in which cell point x is
inline int cell_index(const Grid& grid, double x) {
  return std::floor((x - grid.domain.first) / grid.dx) + grid.n_ghost;
}

// Mock grid
inline Grid make_dummy_grid() { return Grid(std::make_pair(0, 1), 102, 1); }

#endif // HYPSYS1D_GRID_HPP
