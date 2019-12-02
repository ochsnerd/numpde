#include <ancse/cfl_condition.hpp>

#include <Eigen/Dense>


/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   double cfl_number) {
  return std::make_shared<FVMCFLCondition>(grid, model, cfl_number);
}

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler,
                   double cfl_number) {
  return std::make_shared<DGCFLCondition>(grid, model, dg_handler, cfl_number);
}
