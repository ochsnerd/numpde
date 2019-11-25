#ifndef HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
#define HYPSYS1D_FVM_RATE_OF_CHANGE_HPP

#include <memory>
#include <Eigen/Dense>

#include <ancse/config.hpp>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

/// Compute the rate of change due to FVM.
/** The semidiscrete approximation of a PDE using FVM is
 *      du_i/dt = - (F_{i+0.5} - F_{i-0.5}) / dx.
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 * @tparam Reconstruction see e.g. `PWConstantReconstruction`.
 */
template <class NumericalFlux, class Reconstruction>
class FVMRateOfChange : public RateOfChange {
public:
  using Vector = typename Model::Vector;
  using Matrix = typename RateOfChange::Matrix;

  FVMRateOfChange(const Grid &grid,
                  const std::shared_ptr<Model> &model,
                  const NumericalFlux &numerical_flux,
                  const Reconstruction &reconstruction)
    : grid(grid),
      model(model),
      numerical_flux(numerical_flux),
      reconstruction(reconstruction) {}

  virtual void operator()(Matrix& dudt,
                          const Matrix& u0) const override {
    int n_cells = grid.n_cells;
    int n_ghost = grid.n_ghost;

    double dx = grid.dx;
    Vector fL(model->get_nvars()), fR(model->get_nvars());
    // very unsure how portable that actually is, but it works for Eigen::DenseBase
    // and for double, so that's cool I guess ¯\_(ツ)_/¯
    fR *= 0.0;

    // quick sanity check
    assert(dudt.rows() == u0.rows());
    assert(dudt.cols() == u0.cols());
    assert(dudt.cols() == n_cells);
    assert(dudt.rows() == model->get_nvars());


    // This uses the example interface
    reconstruction.set(u0);

    for (int i = n_ghost - 1; i < n_cells - n_ghost; ++i) {
      // This is for when we use real reconstructions
      //auto [uL, uR] = reconstruction(u0, i);

      // This uses the example interface
      auto [uL, uR] = reconstruction(i);

      fL = fR;
      fR = numerical_flux(uL, uR);

      dudt.col(i) = (fL - fR) / dx;
    }
  }

private:
  Grid grid;
  std::shared_ptr<Model> model;
  NumericalFlux numerical_flux;
  Reconstruction reconstruction;
};

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const nlohmann::json &config,
                        const Grid &grid,
                        const std::shared_ptr<Model> &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
