#pragma once

// "Model" stuff
// Also fml, I only just now realized all this stuff is implemented already in euler.hpp
double pressure(const EulerState& u) {
  double gamma = 1.4;
  double p = (u[3] / u[0] - .5 * (u[1] * u[1] + u[2] * u[2]) / (u[0] * u[0])) * (gamma - 1);
  assert(p >= 0);
  return p;
}

double soundspeed(const EulerState& u) {
  double gamma = 1.4;
  double c = std::sqrt(gamma * pressure(u) / u[0]);
  assert(c >= 0);
  return c;
}

EulerState euler_flux(const EulerState& u) {
  auto p = pressure(u);
  return EulerState(
    u[1],
    u[1] * u[1] / u[0] + p,
    u[1] * u[2] / u[0],
    u[1] / u[0] * (u[3] + p));
}

// ASSUMING WE'RE JUST INTERESTED IN WHATS GOING ON THE X DIRECTION SINCE WE
// ALIGN THE CELL BOUNDARY PARALLEL TO THE Y AXIS
double min_euler_eigenvalue(const EulerState& u) {
  return u[1] / u[0] - soundspeed(u);
}

double max_euler_eigenvalue(const EulerState& u) {
  return u[1] / u[0] + soundspeed(u);
}

// intermediate state
EulerState uStar(const EulerState& u, double s, double sStar) {
  // Toro 10.39
  double drho = u[0] * (s - u[1] / u[0]) / (s - sStar);
  return drho * EulerState(
    1,
    sStar,
    u[2] / u[0],
    u[3] / u[0] + (sStar - u[1] / u[0]) * (sStar + pressure(u) / (u[0] * (s - u[1] / u[0]))));
}

EulerState intermediate_flux(double s, const EulerState& uStar, const EulerState& u) {
  return euler_flux(u) + s * (uStar - u);
}

// Speeds (not actually Einfeldt-Batten, sry)
double sL(const EulerState& uL, const EulerState& uR) {
  return std::min(min_euler_eigenvalue(uL), min_euler_eigenvalue(uR));
}

double sR(const EulerState& uL, const EulerState& uR) {
  return std::max(max_euler_eigenvalue(uL), max_euler_eigenvalue(uR));
}

double sStar(const EulerState& uL, const EulerState& uR, double sL, double sR) {
  // Toro 10.37
  auto rhoL = uL[0];
  auto rhoR = uR[0];
  auto pL = pressure(uL);
  auto pR = pressure(uR);
  auto vL = uL[1] / uL[0];
  auto vR = uR[1] / uR[0];

  return (pR - pL
          + rhoL * vL * (sL - vL)
          - rhoR * vR * (sR - vR))
         /
         (rhoL * (sL - vL)
          - rhoR * (sR - vR));
}

/// HLLC numerical flux with Einfeldt-Batten wavespeeds.
/** Reference: Batten, Wavespeed Estimates for the HLLC Riemann Solver, 1997
 *  @param euler
 *  @param uL    conserved variables 'inside' of the cell.
 *  @param uR    conserved variables 'outside' of the cell.
 */
EulerState hllc(const EulerState &uL, const EulerState &uR) {
  // implement the HLLC approximate Riemann solver.
  // tip, compute the three wave speeds sL, s_star & sR in
  // a separate function to keep things readable.
  double sL_ = sL(uL, uR);
  if (sL_ > 0) {
    return euler_flux(uL);
  }

  double sR_ = sR(uL, uR);
  if (sR_ < 0) {
    return euler_flux(uR);
  }

  double sStar_ = sStar(uL, uR, sL_, sR_);
  if (sStar_ > 0) {
    return intermediate_flux(sL_, uStar(uL, sL_, sStar_), uL);
  }

  return intermediate_flux(sR_, uStar(uR, sR_, sStar_), uR);
}


// EulerState rusanov(const EulerState& uL, const EulerState& uR, const Eigen::Vector2D& n) {
//   double speed = std::max(max_euler_eigenvalue(uL, n), max_euler_eigenvalue(uR, n));

//   return .5 * (euler_flux(uR) + euler_flux(uL)) - .5 * speed * (uR -uL);
// }

// double max_euler_eigenvalue(const EulerState& u, const Eigen::Vector2D& n) {
//   double gamma = 1.4;
//   double p = (u[3] / u[0] - .5 * (u[1] * u[1] + u[2] * u[2])) * (gamma - 1);
//   double c = std::sqrt(gamma * p / u[0]);
//   return std::abs(n.dot(Eigen::Vector2D(u[1] / u[0], u[2] / u[0]))) + c;
// }

// EulerState euler_flux(const EulerState& u) {
  
// }
