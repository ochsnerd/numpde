#include <ancse/dg_rate_of_change.hpp>

#include <Eigen/Dense>
#include <ancse/config.hpp>
#include <ancse/polynomial_basis.hpp>
#include <ancse/dg_handler.hpp>
#include <ancse/numerical_flux.hpp>
#include <fmt/format.h>


/// DG numerical flux term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_numerical_flux (Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
  int n_cells = grid.n_cells;
  int n_ghost = grid.n_ghost;
  int n_vars = model->get_nvars();

  double dx = grid.dx;
  Eigen::VectorXd fL(n_vars), fR(n_vars);

  fR *= 0;

  // quick sanity check
  assert(dudt.rows() == u0.rows());
  assert(dudt.cols() == u0.cols());
  assert(dudt.cols() == n_cells);
  assert(dudt.rows() == n_vars);

  for (int i = n_ghost - 1; i < n_cells - n_ghost; ++i) {
    fL = fR;

    // I know this is less than ideal
    auto uL = dg_handler.build_sol(u0.col(i), cell_center(grid, i) - 0.49 * grid.dx);
    auto uR = dg_handler.build_sol(u0.col(i), cell_center(grid, i) + 0.49 * grid.dx);
    fR = numerical_flux(uL, uR);

    dudt.col(i) = (fL - fR) / dx;
  }
}

/// DG volume integral term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_volume_integral(Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
    // implement the loop for DG volume integral.
}


#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)          \
    if (config["flux"] == (token)) {                            \
        return std::make_shared< DGRateOfChange<FluxType> >(    \
            grid, model, flux, poly_basis, dg_handler);                     \
    }


std::shared_ptr<RateOfChange> make_dg_rate_of_change(
    const nlohmann::json &config,
    const Grid &grid,
    const std::shared_ptr<Model> &model,
    const PolynomialBasis &poly_basis,
    const DGHandler &dg_handler,
    const std::shared_ptr<SimulationTime> &simulation_time)
{
    // Register the other numerical fluxes.
  REGISTER_NUMERICAL_FLUX("Rusanov_flux", RusanovFlux, RusanovFlux(model))

    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}",
                    std::string(config["flux"])));
}

#undef REGISTER_NUMERICAL_FLUX
