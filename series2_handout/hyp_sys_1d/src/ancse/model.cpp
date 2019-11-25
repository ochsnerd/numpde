#include <ancse/model.hpp>

#include <iostream>
#include <fmt/format.h>


///------------------///
/// Euler equations  ///
///------------------///

#define REGISTER_MODEL(token, ModelType)      \
    if (config["model"] == (token)) {         \
        return std::make_shared<ModelType>(); \
    }

std::shared_ptr<Model> make_model (const nlohmann::json &config)
{
    REGISTER_MODEL("burgers", Burgers)

    REGISTER_MODEL("euler", Euler)

    throw std::runtime_error(
        fmt::format("Unknown model. {}", std::string(config["flux"])));
}

#undef REGISTER_MODEL
