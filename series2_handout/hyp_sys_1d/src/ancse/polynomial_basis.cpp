#include "../../include/ancse/polynomial_basis.hpp"

#include <cmath>


// implementing the derivative of the Legendre polynomials
// for general q seems a bit overkill (since all explicit
// formulas seem to involve binomial coefficents and I don't
// feel like dealing with that).
// Also I'm a dumbass and wasn't working with the most up-to-date
// exercise sheet.
double legendre_deriv(int k, double xi) {
  assert(k >= 0 && k <= 4);
  switch (k) {
  case 0 : return 0;
  case 1 : return 1;
  case 2 : return 3 * xi;
  case 3 : return .5 * (15 * xi * xi - 3);
  case 4 : return 2.5 * (7 * xi * xi * xi - 3 * xi);
  }

  // We should never end up here (as long as assertions are on)
  throw std::invalid_argument("0 <= k <= 4");
  return 1;
}

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
  assert(0 <= xi && xi <= 1);

  Eigen::VectorXd dphi(p + 1);

  for (int k = 0; k < p + 1; ++k) {
    dphi(k) = std::sqrt(2 * k + 1) * legendre_deriv(k, 2 * xi - 1);
  }

  return scaling_factor * dphi;
}
