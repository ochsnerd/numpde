#include <Eigen/Dense>

#include <ancse/cfl_condition.hpp>
#include <ancse/config.hpp>
#include <ancse/fvm_rate_of_change.hpp>
#include <ancse/snapshot_writer.hpp>
#include <ancse/time_loop.hpp>

template<class F>
Eigen::VectorXd ic(const F &f, const Grid &grid) {
    Eigen::VectorXd u0(grid.n_cells);
    for(int i = 0; i < grid.n_cells; ++i) {
        u0[i] = f(cell_center(grid, i));
    }

    return u0;
}

TimeLoop make_fvm(const Grid &grid) {
    auto config = get_global_config();
    double t_end = config["t_end"];
    double cfl_number = config["cfl_number"];

    auto n_ghost = grid.n_ghost;
    auto n_cells = grid.n_cells;

    auto model = Model{};

    auto simulation_time = std::make_shared<SimulationTime>(t_end);
    auto fvm_rate_of_change = make_fvm_rate_of_change(grid, model, simulation_time);
    auto boundary_condition = make_boundary_condition(n_ghost, config["boundary_condition"]);
    auto time_integrator = make_runge_kutta(fvm_rate_of_change, boundary_condition, n_cells);
    auto cfl_condition = make_cfl_condition(grid, model, cfl_number);
    auto snapshot_writer = std::make_shared<JSONSnapshotWriter>(grid, simulation_time, std::string(config["output"]));

    return TimeLoop(simulation_time, time_integrator, cfl_condition, snapshot_writer);
}

void smooth_sine_test() {
    auto config = get_global_config();

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({0.0, 1.0}, n_cells, n_ghost);
    auto u0 = ic([](double x) { return std::sin(2.0*M_PI * x);}, grid);

    auto fvm = make_fvm(grid);
    fvm(u0);
}

int main() {
    smooth_sine_test();

    return 0;
}