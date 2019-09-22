#include "writer.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

void apply_boundary_conditions(Eigen::ArrayXXd &u, int k) {
    auto N = u.rows();
    u(0, k) = u(1, k);
    u(N - 1, k) = u(N - 2, k);
}

//----------------upwindFDBegin----------------
/// Uses forward Euler and upwind finite differences to compute u from time 0 to
/// time T
///
/// @param[in] u0 the initial conditions in the physical domain, excluding
/// ghost-points.
/// @param[in] dt the time step size
/// @param[in] T the solution is computed for the interval [0, T]. T is assumed
/// to be a multiple of of dt.
/// @param[in] a the advection velocity
/// @param[in] domain left & right limit of the domain
///
/// @return returns the solution 'u' at every time-step and the corresponding
/// time-steps. The solution `u` includes the ghost-points.
std::pair<Eigen::ArrayXXd, Eigen::VectorXd>
upwindFD(const Eigen::VectorXd &u0,
         double dt,
         double T,
         const std::function<double(double)> &a,
         const std::pair<double, double> &domain) {

    auto N = u0.size();
    auto nsteps = int(round(T / dt));

    auto u = Eigen::ArrayXXd(N + 2, nsteps + 1);
    auto time = Eigen::VectorXd(nsteps + 1);

    auto [xL, xR] = domain;
    double dx = (xR - xL) / (N - 1.0);

    // initial conditions
    u.block(1, 0, N, 1) = u0;

    // outflow boundary conditions
    apply_boundary_conditions(u, 0);

    /* Main loop */
    for (int t = 0; t < nsteps; ++t) {
      time(t) = t * dt;
      for (int i = 1; i < N + 1; ++i) {
        auto a_ = a(xL + (i - 1) * dx);

        auto a_plus = std::max(a_, 0.);
        auto a_minus = std::min(a_, 0.);

        auto u_plus = (u(i + 1, t) - u(i, t)) / dx;
        auto u_minus = (u(i, t) - u(i - 1, t)) / dx;

        u(i, t + 1) = u(i, t) - dt * (a_plus * u_minus + a_minus * u_plus);
      }
      apply_boundary_conditions(u, t + 1);
    }
    time(nsteps) = T;

    return {std::move(u), std::move(time)};
}
//----------------upwindFDEnd----------------

//----------------centeredFDBegin----------------
/// Uses forward Euler and centered finite differences to compute u from time 0
/// to time T
///
/// @param[in] u0 the initial conditions, as column vector
/// @param[in] dt the time step size
/// @param[in] T the solution is computed for the interval [0, T]. T is assumed
/// to be a multiple of of dt.
/// @param[in] a the advection velocity
/// @param[in] domain left & right limit of the domain
///
/// @return returns the solution 'u' at every time-step and the corresponding
/// time-steps. The solution `u` includes the ghost-points.
std::pair<Eigen::ArrayXXd, Eigen::VectorXd>
centeredFD(const Eigen::VectorXd &u0,
           double dt,
           double T,
           const std::function<double(double)> &a,
           const std::pair<double, double> &domain) {

    auto N = u0.size();
    auto nsteps = int(round(T / dt));
    auto u = Eigen::ArrayXXd(N + 2, nsteps + 1);
    auto time = Eigen::VectorXd(nsteps + 1);

    auto [xL, xR] = domain;
    double dx = (xR - xL) / (N - 1.0);

    // initial conditions
    u.block(1, 0, N, 1) = u0;

    // outflow boundary conditions
    apply_boundary_conditions(u, 0);

    /* Main loop */
    for (int t = 0; t < nsteps; ++t) {
      time(t) = t * dt;
      for (int i = 1; i < N + 1; ++i) {
        u(i, t + 1) = u(i + 1, t) - u(i - 1, t);
        u(i, t + 1) *= (-.5) * dt * a(xL + (i - 1) * dx) / dx;
        u(i, t + 1) += u(i, t);
      }
      apply_boundary_conditions(u, t + 1);
    }
    time(nsteps) = T;

    return {std::move(u), std::move(time)};
}
//----------------centeredFDEnd----------------

/* Initial condition: rectangle */
double ic(double x) {
    if (x < 0.25 || x > 0.75)
        return 0.0;
    else
        return 2.0;
}

int main() {
    double T = 2.0;
    double dt = 0.002; // Change this for timestep comparison
    int N = 101;

    double xL = 0.0;
    double xR = 5.0;
    auto domain = std::pair<double, double>{xL, xR};

    auto a = [](double x) { return std::sin(2.0 * M_PI * x); };
    //auto a = [](double x) { return 1; };

    Eigen::VectorXd u0(N);
    double h = (xR - xL) / (N - 1.0);
    /* Initialize u0 */
    for (int i = 0; i < u0.size(); i++) {
        u0[i] = ic(xL + h * i);
    }

    const auto &[u_upwind, time_upwind] = upwindFD(u0, dt, T, a, domain);
    writeToFile("time_upwind.txt", time_upwind);
    writeMatrixToFile("u_upwind.txt", u_upwind);

    const auto &[u_centered, time_centered] = centeredFD(u0, dt, T, a, domain);
    writeToFile("time_centered.txt", time_centered);
    writeMatrixToFile("u_centered.txt", u_centered);
}
