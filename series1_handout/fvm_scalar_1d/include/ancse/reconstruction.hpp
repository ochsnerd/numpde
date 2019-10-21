#ifndef FVMSCALAR1D_RECONSTRUCTION_HPP
#define FVMSCALAR1D_RECONSTRUCTION_HPP

#include <Eigen/Dense>
#include <cmath>
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

inline double sign(double a) { return copysign(1.0, a); }

inline double slope(double L, double R, double dx) { return (R - L) / dx; }

inline double minmod(double a, double b) {
  return sign(a) == sign(b) ? sign(a) * std::min(std::abs(a),std::abs(b)) : 0;
}

inline double maxmod(double a, double b){
  return sign(a) == sign(b) ? sign(a) * std::max(std::abs(a),std::abs(b)) : 0;
}

class PWConstantReconstruction {
  public:
    /// Compute the left and right trace at the interface i + 1/2.
    /** Note: This API is agnostic to the number of cell-averages required
     *        by the method. Therefore, reconstructions with different stencil
     *        sizes can implement this API; and this call can be used in parts
     *        of the code that do not need to know about the details of the
     *        reconstruction.
     */
    inline std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                                int i) const {
        return (*this)(u[i], u[i + 1]);
    }

    /// Compute the left and right trace at the interface.
    /** Piecewise constant reconstruction of the left and right trace only
     *  requires the cell-average to the left and right of the interface.
     *
     *  Note: Compared to the other overload this reduces the assumption on
     *        how the cell-averages are stored. This is useful when testing and
     *        generally makes the function useful in more situations.
     */
    inline std::pair<double, double> operator()(double ua, double ub) const {
        return {ua, ub};
    }
};

template <class SlopeLimiter>
class SecondOrderReconstruction {
public:
  SecondOrderReconstruction(const Grid& grid, SlopeLimiter sigma) : dx_{grid.dx}, sigma_{sigma} {}

  std::pair<double, double> operator() (const Eigen::VectorXd& u, int i) const {
    // unsure if correct, what is actually stored in u?
    // u_i or u_{i + 1/2}?
    return {u[i] + .5 * dx_ * sigma_(u, dx_, i), u[i + 1] - .5 * dx_ * sigma_(u, dx_, i + 1)};
  }

private:
  double dx_;
  SlopeLimiter sigma_;
};

inline double minmod_limiter(const Eigen::VectorXd& u, double dx, int i) {
  return minmod(slope(u[i], u[i+1], dx), slope(u[i-1], u[i], dx));
}

inline double superbee_limiter(const Eigen::VectorXd& u, double dx, int i) {
  auto slope_l = slope(u[i-1], u[i], dx);
  auto slope_r = slope(u[i], u[i+1], dx);

  auto sigma_l = minmod(2*slope_l, slope_r);
  auto sigma_r = minmod(slope_l, 2*slope_r);

  return maxmod(sigma_l, sigma_r);
}

#endif // FVMSCALAR1D_RATE_OF_CHANGE_HPP
