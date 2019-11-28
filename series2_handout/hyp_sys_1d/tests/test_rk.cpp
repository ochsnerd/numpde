// So this is basically copy-pasted from the solutions to series 1
#include <ancse/runge_kutta.hpp>

#include <ancse/boundary_condition.hpp>
#include <fmt/format.h>
#include <gtest/gtest.h>

/// Compute the convergence rate.
/** This assumes that the controling parameter, e.g. number of steps, number of
 * cells, mesh width, step-size is halved/double between u_coarse and u_fine.
 */
double convergence_rate(double err_coarse, double err_fine) {
    return (std::log(err_coarse) - std::log(err_fine)) / std::log(2);
}

double almost_equal(double a, double b, double atol) {
    return std::abs(a - b) < atol;
}

/// Right hand side of the ODE: dudt = -u.
/** This simple RateOfChange can be used to check if SSP2
 *  converges at the right order.
 */
class ExpODE : public RateOfChange {
  public:
    virtual void operator()(Eigen::MatrixXd &dudt,
                            const Eigen::MatrixXd &u0) const override {
        dudt = -u0;
    }

  virtual void to_cons(Eigen::MatrixXd& u) const override {}
};

double
solve_ode(const TimeIntegrator &rk, double ic, double t_end, int n_steps) {
  Eigen::MatrixXd u0(1,1);
  u0(0, 0) = ic;

  double dt = t_end / n_steps;

  Eigen::MatrixXd u1(1,1);
  for (int i = 0; i < n_steps; ++i) {
    rk(u1, u0, dt);

    u1.swap(u0);
  }

  return u0(0,0);
}

// This an example of how to mocks can be useful.
//
// If we run SSP2 with the semi-discrete formulation of FVM as a right-hand side
// and observe that the rate is not 2, we've not learnt much. Just that
// somewhere in our code there is a bug. This is also useful, and would make a
// *very* important test; but should be implemented elsewhere.
//
// Here we want to only test SSP2 (or RK schemes in general). Since it is
// implemented against a generic RHS (called RateOfChange), we can also solve
// the ODE:
//     d/dt u = -u
// We know the exact solution and can compute convergence rates. Furthermore,
// the problem gets isolated to code in this file and code strictly related to
// SSP2, nothing else.
//
// ExpODE is what's referred to as a 'mock'.
//
// NOTE: Asserting theoretical convergence rates is very powerful if the rate is
//       at least 2. Many small mistakes have the unfortunate property of still
//       being first order accurate, e.g. off-by-one indexing mistakes.
//
// NOTE: This test accept a reference to the abstract base class. You can run
//       this for any RK scheme.
void check_rk_on_exp_ode(const TimeIntegrator &rk, double expected_rate) {
    int n_steps = 20;
    double t_end = 1.0;

    double u0 = 2.0;

    double u_coarse = solve_ode(rk, u0, t_end, n_steps);
    double u_fine = solve_ode(rk, u0, t_end, 2 * n_steps);
    double u_exact = 2.0 * std::exp(-t_end);

    double err_coarse = std::abs(u_coarse - u_exact);
    double err_fine = std::abs(u_fine - u_exact);

    double observed_rate = convergence_rate(err_coarse, err_fine);

    double atol = 0.1;
    ASSERT_TRUE(almost_equal(observed_rate, expected_rate, atol))
        << fmt::format("rate = %.2f, err_coarse = %.3e, err_fine = %.3e",
                       observed_rate,
                       err_coarse,
                       err_fine);
}

TEST(TimeIntegrator, ExpODE) {
    int n_cells = 1;
    int n_ghost = 0;

    // We can abuse OutflowBC with n_ghost = 0. Or one could have implemented a
    // NoBC, which just does nothing.
    auto bc = std::make_shared<OutflowBC>(n_ghost);
    auto rhs = std::make_shared<ExpODE>();

    auto forward_euler = ForwardEuler(rhs, bc, nullptr, 1, n_cells);
    check_rk_on_exp_ode(forward_euler, 1.0);

    auto ssp2 = SSPRK2(rhs, bc, nullptr, 1, n_cells);
    check_rk_on_exp_ode(ssp2, 2.0);
}
