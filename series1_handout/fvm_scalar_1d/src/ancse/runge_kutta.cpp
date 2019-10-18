#include <ancse/runge_kutta.hpp>

#include <ancse/config.hpp>
#include <fmt/format.h>

#define REGISTER_RUNGE_KUTTA(token, RKType)                                    \
    if (rk_key == token) {                                                     \
        return std::make_shared<RKType>(                                       \
            rate_of_change, boundary_condition, n_cells);                      \
    }

std::shared_ptr<RungeKutta>
make_runge_kutta(const std::shared_ptr<RateOfChange> &rate_of_change,
                 const std::shared_ptr<BoundaryCondition> &boundary_condition,
                 int n_cells) {
    auto c = get_global_config();
    std::string rk_key = c["time_integrator"];

    REGISTER_RUNGE_KUTTA("forward_euler", ForwardEuler);

    // Register your SSP2 class.

    throw std::runtime_error(
        fmt::format("Unknown time-integrator. [{}]", rk_key));
}
