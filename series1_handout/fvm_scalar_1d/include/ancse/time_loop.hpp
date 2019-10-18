#ifndef FVMSCALAR1D_TIME_LOOP_HPP
#define FVMSCALAR1D_TIME_LOOP_HPP

#include <Eigen/Dense>
#include <fmt/format.h>
#include <fstream>

#include <ancse/cfl_condition.hpp>
#include <ancse/runge_kutta.hpp>
#include <ancse/simulation_time.hpp>
#include <ancse/snapshot_writer.hpp>

/// Advances the solution of the PDE iteratively until the final time.
class TimeLoop {
  public:
    TimeLoop(std::shared_ptr<SimulationTime> simulation_time,
             std::shared_ptr<TimeIntegrator> time_integrator,
             std::shared_ptr<CFLCondition> cfl_condition,
             std::shared_ptr<SnapshotWriter> snapshot_writer);

    void operator()(Eigen::VectorXd u0) const;

    void write_snapshot(const Eigen::VectorXd &u) const;

  private:
    mutable std::shared_ptr<SimulationTime> simulation_time;
    std::shared_ptr<TimeIntegrator> time_integrator;
    std::shared_ptr<CFLCondition> cfl_condition;
    std::shared_ptr<SnapshotWriter> snapshot_writer;
};

#endif // FVMSCALAR1D_TIME_LOOP_HPP
