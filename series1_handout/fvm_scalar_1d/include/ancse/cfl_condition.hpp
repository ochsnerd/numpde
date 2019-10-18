#ifndef FVMSCALAR1D_CFL_CONDITION_HPP
#define FVMSCALAR1D_CFL_CONDITION_HPP

#include <cmath>
#include <limits>
#include <memory>

#include <Eigen/Dense>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>

/// Interface for computing CFL restricted timesteps.
class CFLCondition {
  public:
    virtual ~CFLCondition() = default;

    /// Compute the largest time-step satisfying the CFL condition.
    virtual double operator()(const Eigen::VectorXd &u) const = 0;
};


std::shared_ptr<CFLCondition> make_cfl_condition(const Grid &grid, const Model &model, double cfl_number);

#endif // FVMSCALAR1D_CFL_CONDITION_HPP
