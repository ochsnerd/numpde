#ifndef FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP
#define FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP

#include <Eigen/Dense>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>
#include <memory>

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
    FVMRateOfChange(const Grid &grid,
                    const NumericalFlux &numerical_flux,
                    const Reconstruction &reconstruction)
        : grid(grid),
          numerical_flux(numerical_flux),
          reconstruction(reconstruction) {}

    virtual void operator()(Eigen::VectorXd &dudt,
                            const Eigen::VectorXd &u0) const override {
        // implement the flux loop here.
    }

  private:
    Grid grid;
    NumericalFlux numerical_flux;
    Reconstruction reconstruction;
};

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const Grid &grid,
                        const Model &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP
