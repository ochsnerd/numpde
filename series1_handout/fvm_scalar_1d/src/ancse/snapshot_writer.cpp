// Putting the header as the very first line of code, ensures that your
// header is self-contained. It's impossible for include-ordering to break
// your code.
#include <ancse/snapshot_writer.hpp>

#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <fstream>

JSONSnapshotWriter::JSONSnapshotWriter(
    const Grid &grid,
    std::shared_ptr<SimulationTime> simulation_time,
    std::string basename)
    : grid(grid),
      simulation_time(std::move(simulation_time)),
      basename(std::move(basename)) {}

void JSONSnapshotWriter::operator()(const Eigen::VectorXd &u) const {
    std::string filename = basename + fmt::format("{:04d}.json", k_output);
    k_output += 1;

    auto json = nlohmann::json{};
    for (int i = 0; i < grid.n_cells; ++i) {
        json["data"][i] = u[i];
        json["cell_centers"][i] = cell_center(grid, i);
    }

    json["time"] = simulation_time->t;

    auto file = std::ofstream(filename);
    assert(file.good());  // always check that you can write to the file.

    file << json.dump(2);
}
