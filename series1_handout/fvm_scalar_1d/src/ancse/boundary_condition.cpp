#include <ancse/boundary_condition.hpp>

#include <fmt/format.h>

std::shared_ptr<BoundaryCondition>
make_boundary_condition(int n_ghost, const std::string &bc_key) {

    // register you boundary conditions here.

    throw std::runtime_error(
        fmt::format("Unknown boundary condition. [{}]", bc_key));
}
