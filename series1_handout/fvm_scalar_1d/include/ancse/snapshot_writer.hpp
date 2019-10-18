#ifndef FVMSCALAR1D_SNAPSHOT_WRITER_HPP
#define FVMSCALAR1D_SNAPSHOT_WRITER_HPP

#include <Eigen/Dense>
#include <ancse/grid.hpp>
#include <ancse/simulation_time.hpp>
#include <memory>

/// Interface for outputting snapshots.
class SnapshotWriter {
  public:
    virtual ~SnapshotWriter() = default;

    virtual void operator()(const Eigen::VectorXd &u) const = 0;
};

class JSONSnapshotWriter : public SnapshotWriter {
  public:
    JSONSnapshotWriter(const Grid &grid,
                      std::shared_ptr<SimulationTime> simulation_time,
                       std::string basename);

    virtual void operator()(const Eigen::VectorXd &u) const override;

  private:
    Grid grid;
    std::shared_ptr<SimulationTime> simulation_time;

    std::string basename;
    mutable int k_output = 0;
};

#endif // FVMSCALAR1D_SNAPSHOT_WRITER_HPP
