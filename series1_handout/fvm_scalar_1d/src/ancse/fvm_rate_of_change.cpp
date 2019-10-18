#include <ancse/fvm_rate_of_change.hpp>

#include <ancse/config.hpp>
#include <ancse/numerical_flux.hpp>
#include <ancse/reconstruction.hpp>
#include <fmt/format.h>

#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)                         \
    if (config["flux"] == (token)) {                                           \
        return std::make_shared<FVMRateOfChange<FluxType, Reconstruction>>(    \
            grid, flux, reconstruction);                                       \
    }

template <class Reconstruction>
std::shared_ptr<RateOfChange>
deduce_numerical_flux(const Grid &grid,
                      const Model &model,
                      const std::shared_ptr<SimulationTime> &simulation_time,
                      const Reconstruction &reconstruction) {
    auto config = get_global_config();

    REGISTER_NUMERICAL_FLUX("central_flux", CentralFlux, CentralFlux(model))

    // Register the other numerical fluxes.

    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}", std::string(config["flux"])));
}
#undef REGISTER_NUMERICAL_FLUX

#define REGISTER_RECONSTRUCTION(token, reconstruction)                         \
    if (config["reconstruction"] == token) {                                   \
        return deduce_numerical_flux(                                          \
            grid, model, simulation_time, reconstruction);                     \
    }

std::shared_ptr<RateOfChange> make_fvm_rate_of_change(
    const Grid &grid,
    const Model &model,
    const std::shared_ptr<SimulationTime> &simulation_time) {

    auto config = get_global_config();

    REGISTER_RECONSTRUCTION("o1", PWConstantReconstruction{})

    // Register the other numerical fluxes.

    throw std::runtime_error(fmt::format(
        "Unknown reconstruction. [{}]", std::string(config["reconstruction"])));
}

#undef REGISTER_RECONSTRUCTION
