#include <ancse/dg_handler.hpp>


/// build solution from DG coefficients and the basis
/// pre-evaluated at a certain point
Eigen::VectorXd DGHandler
:: build_sol(const Eigen::VectorXd& u,
             const Eigen::VectorXd& basis) const
{
  // untested!
  assert(u.size() == n_vars * n_coeff);
  assert(basis.size() == n_coeff);

  Eigen::VectorXd u_h(n_vars);

  for (int l = 0; l < n_vars; ++l) {
    u_h(l) = basis.dot(u.segment(l * n_coeff, n_coeff)) / std::sqrt(grid.dx);
  }

  return u_h;
}

/// build solution from DG coefficients at a given reference point
Eigen::VectorXd DGHandler 
:: build_sol(const Eigen::VectorXd& u,
             double xi) const
{
  // untested!
  // well, reference_point is tested and this function doesn't do much
  // more than that...
  double x = reference_point(grid, xi);
  return build_sol(u, poly_basis(x));
}

/// build cell average
Eigen::MatrixXd DGHandler
:: build_cell_avg (const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();

    assert(n_cells == grid.n_cells);
    assert(u.rows() == n_vars * n_coeff);

    Eigen::MatrixXd u0(n_vars, n_cells);

    // I'm sure this could be done cooler with some Eigen
    // fucntionality involving stride
    for (int i = 0; i < n_vars; ++i) {
      for (int j = 0; j < n_cells; ++j) {
        u0(i, j) = u(i * n_coeff, j);
      }
    }

    return u0;
}

/// build split solution uSol_m = u0 + um, uSol_p = u0 - up
/// from DG coefficients
std::tuple <Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
DGHandler :: build_split_sol(const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();
    if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    
    auto u0 = build_cell_avg(u);
    Eigen::MatrixXd um = Eigen::MatrixXd::Zero (n_vars, n_cells);
    Eigen::MatrixXd up = Eigen::MatrixXd::Zero (n_vars, n_cells);


    return {std::move(u0), std::move(um), std::move(up)};
}

/// build DG coefficients from uSol_m = u0 + um, uSol_p = u0 - up
void DGHandler :: compute_limit_coeffs (Eigen::MatrixXd &u,
                                        Eigen::MatrixXd &um,
                                        Eigen::MatrixXd &up) const
{
    if (n_coeff == 1) {
        return;
    }
    else if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    
}
