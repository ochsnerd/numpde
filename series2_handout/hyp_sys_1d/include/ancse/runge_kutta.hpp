#ifndef HYPSYS1D_RUNGE_KUTTA_HPP
#define HYPSYS1D_RUNGE_KUTTA_HPP

#include <Eigen/Dense>
#include <map>
#include <memory>

#include <ancse/includes.hpp>
#include <ancse/config.hpp>
#include <ancse/dg_limiting.hpp>
#include <ancse/boundary_condition.hpp>
#include <ancse/rate_of_change.hpp>

/// Interface advancing the solution of a PDE by one step.
class TimeIntegrator {
public:
  using Matrix = Eigen::MatrixXd;

  virtual ~TimeIntegrator() = default;

  /// Update 'u1' starting from 'u0' using a time-step `dt`.
  virtual void operator()(Matrix& u1,
                          const Matrix& u0,
                          double dt) const = 0;
};

/// Runge-Kutta time integration methods.
/** Note: not all methods are Runge-Kutta methods. Therefore, 'RungeKutta' can
 *        not be the interface. It should be the base for all (explicit)
 *        Runge-Kutta methods.
 */
class RungeKutta : public TimeIntegrator {
  public:
  using Matrix = typename TimeIntegrator::Matrix;

  RungeKutta(std::shared_ptr<RateOfChange> rate_of_change_,
             std::shared_ptr<BoundaryCondition> boundary_condition_,
             std::shared_ptr<Limiting> limiting_)
    : rate_of_change(std::move(rate_of_change_)),
      boundary_condition(std::move(boundary_condition_)),
      limiting(std::move(limiting_))
  {}

protected:
  void post_euler_step(Eigen::MatrixXd &u) const {
    (*boundary_condition)(u);
    if(limiting != nullptr) {
      (*limiting)(u);
    }
    (*boundary_condition)(u);
  }

protected:
  std::shared_ptr<RateOfChange> rate_of_change;
  std::shared_ptr<BoundaryCondition> boundary_condition;
  std::shared_ptr<Limiting> limiting;
};

class ForwardEuler : public RungeKutta {
  private:
    using super = RungeKutta;

  public:
    ForwardEuler(std::shared_ptr<RateOfChange> rate_of_change_,
                 std::shared_ptr<BoundaryCondition> boundary_condition_,
                 std::shared_ptr<Limiting> limiting_,
                 int n_rows,
                 int n_cols)
        : super(std::move(rate_of_change_),
                std::move(boundary_condition_),
                std::move(limiting_)),
          dudt(n_rows, n_cols) {}

    virtual void operator()(Eigen::MatrixXd &u1,
                            const Eigen::MatrixXd &u0,
                            double dt) const override {

        (*rate_of_change)(dudt, u0);
        u1 = u0 + dt * dudt;
        post_euler_step(u1);
    }

  private:
    mutable Eigen::MatrixXd dudt;
};

class SSPRK2 : public RungeKutta {
public:
  using Matrix = typename RungeKutta::Matrix;

  SSPRK2(std::shared_ptr<RateOfChange> rate_of_change,
         std::shared_ptr<BoundaryCondition> boundary_condition,
         std::shared_ptr<Limiting> limiting,
         int n_vars,
         int n_cells)
    : RungeKutta(std::move(rate_of_change),
                 std::move(boundary_condition),
                 std::move(limiting)),
      dudt(n_vars, n_cells) {}

  virtual void operator() (Matrix& u1,
                           const Matrix& u0,
                           double dt) const override {
    // Assuming that the IC already respect the boundary conditions,
    // so we don't need to apply them before the first timestep.
    (*rate_of_change)(dudt, u0);
    u1 = u0 + dt * dudt;
    post_euler_step(u1);

    (*rate_of_change)(dudt, u1);
    u1 = u1 + dt * dudt;
    post_euler_step(u1);

    u1 = 0.5 * (u0 + u1);
  }
  
private:
  mutable Matrix dudt;
};

/// make Runge Kutta for FVM
std::shared_ptr<RungeKutta>
make_runge_kutta(const nlohmann::json &config,
                 const std::shared_ptr<RateOfChange> &rate_of_change,
                 const std::shared_ptr<BoundaryCondition> &boundary_condition,
                 int n_rows,
                 int n_cols);

/// make Runge Kutta for DG
std::shared_ptr<RungeKutta>
make_runge_kutta(const nlohmann::json &config,
                 const std::shared_ptr<RateOfChange> &rate_of_change,
                 const std::shared_ptr<BoundaryCondition> &boundary_condition,
                 const std::shared_ptr<Limiting> &dg_limiting,
                 int n_rows,
                 int n_cols);

#endif // HYPSYS1D_RUNGE_KUTTA_HPP
