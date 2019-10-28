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

class ConcreteCFL : public CFLCondition {
public:
  ConcreteCFL(const Grid& grid, const Model& model, double cfl_number) :
    model_{model},
    dx_{grid.dx},
    cfl_number_{cfl_number} {};

  virtual double operator()(const Eigen::VectorXd& u) const override {
    auto abs_compare = [](double a, double b) {return std::abs(a) < std::abs(b);};
    double max_abs_flux = std::abs(*std::max_element(&(u[0]),
                                                     &(u[0]) + u.size(),
                                                     abs_compare));

    return cfl_number_ * dx_ / max_abs_flux;
  }

private:
  Model model_;
  double dx_;
  double cfl_number_;
};


std::shared_ptr<CFLCondition> make_cfl_condition(const Grid &grid, const Model &model, double cfl_number);

#endif // FVMSCALAR1D_CFL_CONDITION_HPP
