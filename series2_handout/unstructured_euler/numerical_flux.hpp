#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include "gradient.hpp"
#include "hllc.hpp"
#include "mesh.hpp"
#include "slope_limiter.hpp"

// Note: this class will compute the rate of change due to the fluxes.
// Note: the reason we made this a class is that it allows you to allocate
//       buffers, once at the beginning of the simulation. Add these buffers
//       as needed.
class FluxRateOfChange {
  public:
    explicit FluxRateOfChange(int n_cells) {
        // Allocate buffers as needed.
    }

    void operator()(Eigen::MatrixXd &dudt,
                    const Eigen::MatrixXd &u,
                    const Mesh &mesh) const {

        // Compute the rate of change of u.
        // Note: Please use the method `computeFlux` to abstract
        // away the details of computing the flux through a
        // given interface.
        // Note: You can use `assert_valid_flux` to check
        // if what `computeFlux` returns makes any sense.
        // Note: Do not assume `dudt` is filled with zeros.
      assert(u.cols() == mesh.getNumberOfTriangles());

      dudt.resizeLike(u);

      for (int i = 0; i < u.cols(); ++i) {
        dudt.col(i) *= 0;
        for (int k = 0; k < 3; ++k) {
          EulerState flux = computeFlux(u, i, k, mesh);
          assert_valid_flux(mesh, i, k, flux);
          dudt.col(i) += flux;
        }
      }
    }

    void assert_valid_flux(const Mesh &mesh,
                           int i,
                           int k,
                           const EulerState &nF) const {
        // This is mostly for debugging (but also important to check in
        // real simulations!): Make sure our flux contribution is not
        // nan (ie. it is not not a number, ie it is a number)
        if (!euler::isValidFlux(nF)) {
            // clang-format off
            throw std::runtime_error(
                "invalid value detected in numerical flux, " + euler::to_string(nF)
                + "\nat triangle: " + std::to_string(i)
                + "\nedge:        " + std::to_string(k)
                + "\nis_boundary: " + std::to_string(!mesh.isValidNeighbour(i, k)));
            // clang-format on
        }
    }

    /// Compute the flux through the k-th interface of cell i.
    EulerState computeFlux(const Eigen::MatrixXd &U,
                           int i,
                           int k,
                           const Mesh &mesh) const {
        auto boundary_type = mesh.getBoundaryType(i, k);

        if (boundary_type == Mesh::BoundaryType::INTERIOR_EDGE) {
            return computeInteriorFlux(U, i, k, mesh);
        } else {
            if (boundary_type == Mesh::BoundaryType::OUTFLOW_EDGE) {
                return computeOutflowFlux(U, i, k, mesh);
            } else /* boundary_type == Mesh::BoundaryType::WING_EDGE */
            {
                return computeReflectiveFlux(U, i, k, mesh);
            }
        }
    }

    /// Compute the outflow flux through the k-th interface of cell i.
    /** Note: you know that edge k is an outflow edge.
     */
    EulerState computeOutflowFlux(const Eigen::MatrixXd &U,
                                  int i,
                                  int k,
                                  const Mesh &mesh) const {
      // Implement the outflow flux boundary condition.
      auto u_rotated = rotate_vertical(U.col(i), i, k, mesh);
      auto f_rotated = euler::flux(u_rotated);
      auto f = rotate_back(f_rotated, i, k, mesh);
      return f * mesh.getEdgeLength(i, k);
    }

    /// Compute the reflective boundary flux through the k-th edge of cell i.
    /** Note: you know that edge k is a reflective/wall boundary edge.
     */
    EulerState computeReflectiveFlux(const Eigen::MatrixXd &U,
                                     int i,
                                     int k,
                                     const Mesh &mesh) const {

      // Implement the reflective flux boundary condition.
      auto ui_reconstructed = reconstruction(U.col(i), i, k, mesh);

      auto n = mesh.getUnitNormal(i, k).normalized();
      auto t = Eigen::Vector2d(n[1], -n[0]);
      EulerState uStar = ui_reconstructed;
      // Pretty sure this expression could be simplified considerably because we
      // rotate uStar right after. 
      uStar.segment(1,2) = ui_reconstructed[0] * (-ui_reconstructed.segment(1,2).dot(n) * n + ui_reconstructed.segment(1,2).dot(t) * t);
        

      auto uL_rotated = rotate_vertical(ui_reconstructed, i, k, mesh);
      auto uR_rotated = -rotate_vertical(uStar, i, k, mesh);

      auto f_rotated = hllc(uL_rotated, uR_rotated);
      auto f = rotate_back(f_rotated, i, k, mesh);
      return f * mesh.getEdgeLength(i, k);
    }

    /// Compute the flux through the k-th interface of cell i.
    /** Note: This edge is an interior edge, therefore approximate the flux
     * through this edge with the appropriate FVM formulas.
     */
    EulerState computeInteriorFlux(const Eigen::MatrixXd &U,
                                   int i,
                                   int k,
                                   const Mesh &mesh) const {
      // Reconstruct the trace values of U and compute
      // the numerical flux through the k-th interface of
      // cell i.
      int j = mesh.getNeighbour(i, k);
      assert(mesh.isValidNeighbour(j, k));

      auto ui_reconstructed = reconstruction(U.col(i), i, k, mesh);
      auto uj_reconstructed = reconstruction(U.col(j), j, k, mesh);

      auto uL_rotated = rotate_vertical(ui_reconstructed, i, k, mesh);
      auto uR_rotated = rotate_vertical(uj_reconstructed, j, k, mesh);

      assert(-uR_rotated == rotate_vertical(uj_reconstructed, i, k, mesh));

      auto f_rotated = hllc(uL_rotated, uR_rotated);
      auto f = rotate_back(f_rotated, i, k, mesh);
      return f * mesh.getEdgeLength(i, k);
    }


  private:
  EulerState reconstruction(const EulerState& u, int i, int k, const Mesh& mesh) const {
    // no reconstruction for now
    return u;
  }

  EulerState outflux_through(const Eigen::MatrixXd& u, int i, int k, const Mesh& mesh) const {
    auto u_rotated = rotate_vertical(u.col(i), i, k, mesh);
    auto f_rotated = euler::flux(u_rotated);
    auto f = rotate_back(f_rotated, i, k, mesh);
    return f / mesh.getEdgeLength(i, k);
  }

  EulerState rotate_vertical(EulerState u, int i, int k, const Mesh& mesh) const {
    u.segment(1,2) = rotation(mesh, i, k) * u.segment(1,2);
    return u;
  }

  EulerState rotate_back(EulerState f, int i, int k, const Mesh& mesh) const {
    f.segment(1,2) = invert_rotation(mesh, i, k) * f.segment(1,2);
    return f;
  }

  Eigen::Matrix2d rotation(const Mesh& mesh, int i, int k) const {
    // The method name and documentation really don't agree here...
    auto n = mesh.getUnitNormal(i, k).normalized();
    auto t = Eigen::Vector2d(n[1], -n[0]);

    Eigen::Matrix2d R;
    R.row(0) = n;
    R.row(1) = t;

    return R / sqrt(2);
  }

  Eigen::Matrix2d invert_rotation(const Mesh& mesh, int i, int k) const {
    return rotation(mesh, i, k).transpose();
  }
    // add any member variables you might need here.
};
