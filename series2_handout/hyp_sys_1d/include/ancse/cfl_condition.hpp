#ifndef HYPSYS1D_CFL_CONDITION_HPP
#define HYPSYS1D_CFL_CONDITION_HPP

#include <cmath>
#include <limits>
#include <memory>

#include <Eigen/Dense>

#include <ancse/includes.hpp>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/dg_handler.hpp>

/// Interface for computing CFL restricted timesteps.
class CFLCondition {
  public:
    virtual ~CFLCondition() = default;

    /// Compute the largest time-step satisfying the CFL condition.
    virtual double operator()(const Eigen::MatrixXd &u) const = 0;
};


class FVMCFLCondition : public CFLCondition {
public:
  // Actually already the interface should use this, but then we have
  // a unwanted dependency (which we already have now - CFLCondition is
  // forced to use the same Matrix type as Model)
  using Matrix = Model::Matrix;

  FVMCFLCondition(const Grid& grid,
                       const std::shared_ptr<Model>& model,
                       double cfl_number)
    : grid_{grid},
      model_{model},
      cfl_number_{cfl_number}
  {}

  virtual double operator()(const Matrix& u) const override {
    auto n_cells = grid_.n_cells;
    auto n_ghost = grid_.n_ghost;

    double a_max = 0.0;
    for (int i = grid_.n_ghost; i < n_cells - n_ghost; ++i) {
      a_max = std::max(a_max, model_->max_eigenvalue(u.col(i)));
    }

    return cfl_number_ * grid_.dx / a_max;
  }

private:
  Grid grid_;
  std::shared_ptr<Model> model_;
  double cfl_number_;
};

class DGCFLCondition : public CFLCondition {
public:
  using Matrix = Model::Matrix;

  DGCFLCondition(const Grid& grid,
                 const std::shared_ptr<Model>& model,
                 const DGHandler& dg_handler,
                 double cfl_number)
    : grid_{grid},
      model_{model},
      dg_handler_{dg_handler},
      cfl_number_{cfl_number}
  {}

  virtual double operator()(const Matrix& u) const override {
    auto cell_avgs = dg_handler_.build_cell_avg(u); 
    auto n_cells = grid_.n_cells;
    auto n_ghost = grid_.n_ghost;

    double a_max = 0.0;
    for (int i = grid_.n_ghost; i < n_cells - n_ghost; ++i) {
      a_max = std::max(a_max, model_->max_eigenvalue(cell_avgs.col(i)));
    }

    return cfl_number_ * grid_.dx / a_max;
  }

private:
  Grid grid_;
  std::shared_ptr<Model> model_;
  DGHandler dg_handler_;
  double cfl_number_;
};

/// make CFL condition for FVM
 std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   double cfl_number);

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler,
                   double cfl_number);

#endif // HYPSYS1D_CFL_CONDITION_HPP
