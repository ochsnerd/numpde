#include <ancse/time_loop.hpp>

#include <iostream>

TimeLoop::TimeLoop(std::shared_ptr<SimulationTime> simulation_time,
                   std::shared_ptr<TimeIntegrator> time_integrator,
                   std::shared_ptr<CFLCondition> cfl_condition,
                   std::shared_ptr<SnapshotWriter> snapshot_writer)

    : simulation_time(std::move(simulation_time)),
      time_integrator(std::move(time_integrator)),
      cfl_condition(std::move(cfl_condition)),
      snapshot_writer(std::move(snapshot_writer)) {}

void TimeLoop::operator()(Eigen::VectorXd u0) const {
    Eigen::VectorXd u1(u0.size());

    simulation_time->dt = (*cfl_condition)(u0);

    while (!is_finished(*simulation_time)) {
        (*time_integrator)(u1, u0, simulation_time->dt);

        simulation_time->advance();
        simulation_time->dt = (*cfl_condition)(u1);

        write_snapshot(u1);

        u0.swap(u1);

        std::cout << fmt::format("{: 4d}: {} {}\n",
                                 simulation_time->k,
                                 simulation_time->t,
                                 simulation_time->dt);
    }
}

void TimeLoop::write_snapshot(const Eigen::VectorXd &u) const {
    (*snapshot_writer)(u);
}
