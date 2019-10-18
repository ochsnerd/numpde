#include <ancse/config.hpp>
#include <fstream>

nlohmann::json get_global_config() {
    static nlohmann::json config = []() {
        // The path is relative to the binary. This assumes you're building in
        // in a subfolder and therefore, the config file is located in the
        // parent.
        std::ifstream file("../config.json");
        assert(file.good());

        nlohmann::json config;
        file >> config;
        return config;
    }();

    return config;
}