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


#endif // FVMSCALAR1D_RATE_OF_CHANGE_HPP
