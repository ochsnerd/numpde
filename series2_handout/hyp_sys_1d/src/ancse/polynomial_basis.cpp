#include "../../include/ancse/polynomial_basis.hpp"

#include <cmath>


/// Computes the Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: operator() (double xi) const
{
  assert(0 <= xi && xi <= 1);

  Eigen::VectorXd phi(p + 1);

  for (int k = 0; k < p + 1; ++k) {
    phi(k) = std::sqrt(2 * k + 1) * std::legendre(k, 2 * xi - 1);
  }

  return scaling_factor * phi;
}

/// Computes the derivative of Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: deriv (double xi) const
{
    return Eigen::VectorXd::Zero(p+1);
}
