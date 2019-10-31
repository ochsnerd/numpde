#include <ancse/cfl_condition.hpp>


std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const Model &model, double cfl_number) {
  return std::make_shared<ConcreteCFL>(grid, model, cfl_number);
}
