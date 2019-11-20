#ifndef HYPSYS1D_RECONSTRUCTION_HPP
#define HYPSYS1D_RECONSTRUCTION_HPP

#include <Eigen/Dense>
#include <cmath>
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/limiters.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

template <typename T>
inline T sign(T a) { return copysign(1.0, a); }

template <typename T>
inline double slope(T L, T R, T dx) { return double(R - L) / dx; }

template <typename T>
inline T minmod(T a, T b) {
  return sign(a) == sign(b) ? sign(a) * std::min(std::abs(a),std::abs(b)) : 0;
}

class PWConstantReconstruction {
  public:
    void set(const Eigen::MatrixXd &u) const {
        up.resize(u.rows(),u.cols());
        up = u;
    }

    /// Compute the left and right trace at the interface i + 1/2.
    /** Note: This API is agnostic to the number of cell-averages required
     *        by the method. Therefore, reconstructions with different stencil
     *        sizes can implement this API; and this call can be used in parts
     *        of the code that do not need to know about the details of the
     *        reconstruction.
     */
    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    operator()(int i) const {
        return (*this)(up.col(i), up.col(i+1));
    }

    /// Compute the left and right trace at the interface.
    /** Piecewise constant reconstruction of the left and right trace only
     *  requires the cell-average to the left and right of the interface.
     *
     *  Note: Compared to the other overload this reduces the assumption on
     *        how the cell-averages are stored. This is useful when testing and
     *        generally makes the function useful in more situations.
     */
    inline
    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    operator()(Eigen::VectorXd ua, Eigen::VectorXd ub) const {
        return {std::move(ua), std::move(ub)};
    }

private:
  mutable Eigen::MatrixXd up;
};


template <class ScalarLimiter>
class VectorLimiter {
public:
  VectorLimiter(ScalarLimiter sigma) : sigma_(sigma) {}

  Eigen::VectorXd operator() (const Eigen::MatrixXd& u, double dx, int i) const {
    auto slopes = Eigen::VectorXd(u.rows());
    for (int j = 0; j < u.rows(); ++j) {
      slopes(j) = sigma_(u.row(j), dx, i);
    }

    return slopes;
  }

private:
  ScalarLimiter sigma_;
};


template <class ScalarLimiter>
class AffineReconstruction {
public:
  using Vector = Eigen::VectorXd;

  AffineReconstruction(const Grid& grid, ScalarLimiter sigma) : dx_{grid.dx}, sigma_{sigma} {}

  std::pair<Vector, Vector> operator() (const Eigen::MatrixXd& u, int i) const {
    return {u.col(i) + .5 * dx_ * sigma_(u, dx_, i),
            u.col(i + 1) - .5 * dx_ * sigma_(u, dx_, i + 1)};
  }

private:
  double dx_;
  VectorLimiter<ScalarLimiter> sigma_;
};


inline double minmod_limiter(const Eigen::VectorXd& u, double dx, int i) {
  return minmod(slope(u[i], u[i+1], dx), slope(u[i-1], u[i], dx));
}

#endif // HYPSYS1D_RATE_OF_CHANGE_HPP
