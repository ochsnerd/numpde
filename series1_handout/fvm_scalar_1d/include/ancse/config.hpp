#ifndef FVMSCALAR1D_CONFIG_HPP
#define FVMSCALAR1D_CONFIG_HPP

#include <nlohmann/json.hpp>

/// Returns the config stored at `config.json`.
// Note: this is easily misused as essentially a global variable.
nlohmann::json get_global_config();

#endif // FVMSCALAR1D_CONFIG_HPP
