#ifndef FVMSCALAR1D_NUMERICAL_FLUX_HPP
#define FVMSCALAR1D_NUMERICAL_FLUX_HPP

#include <memory>
#include <cmath>
#include <algorithm>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model. It is also unconditionally a
 * bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * (fL + fR);
    }
  double operator() (std::pair<double,double> u) const {
    return (*this)(u.first, u.second);
  }

  private:
    Model model;
};

template<class ConcreteFlux>
class NumericalFlux {
  // Interface for numerical fluxes

  // Using CRTP to get compile time polymorphism
public:
  virtual ~NumericalFlux() {};

  virtual double compute_flux(double, double) const = 0;

  double operator() (double uL, double uR) const {
    return static_cast<ConcreteFlux const*>(this)->compute_flux(uL, uR);
  }

  double operator() (std::pair<double,double> u) const {
    return static_cast<ConcreteFlux const*>(this)->compute_flux(u.first, u.second);
  }
};

class RusanovFlux : public NumericalFlux<RusanovFlux> {
  // Mishra Hyperbolic PDES 4.2.4
public:
  explicit RusanovFlux(const Model& model) : model_{model} {}

  double compute_flux(double uL, double uR) const override
  {
    auto speed = std::max(std::abs(model_.max_eigenvalue(uL)),
                          std::abs(model_.max_eigenvalue(uR)));

    return CentralFlux(model_)(uL, uR) - 0.5 * (uR - uL) * speed;
  }

private:
  Model model_;
};


class LaxFriedrichsFlux : public NumericalFlux<LaxFriedrichsFlux> {
  // LN 4.2.3
public:
  LaxFriedrichsFlux(const Model& model,
                    const Grid& grid,
                    const std::shared_ptr<SimulationTime>& time) :  // shared_ptr (┛ಠ_ಠ)┛彡┻━┻
    model_{model},
    speed_{grid.dx / time->dt} // shared_ptr (┛ಠ_ಠ)┛彡┻━┻
  {}

  virtual double compute_flux (double uL, double uR) const override {
    return CentralFlux(model_)(uL, uR) - 0.5 * speed_ * (uR -uL);
  }

private:
  Model model_;
  double speed_;
};

class RoeFlux : public NumericalFlux<RoeFlux> {
  // Mishra Hyperbolic PDES 4.2.1
public:
  explicit RoeFlux(const Model& model) : model_{model} {}

  virtual double compute_flux(double uL, double uR) const override {
    double linearized_flux;
    if (std::abs(uL - uR) > 1e-12) {
      linearized_flux = (model_.flux(uR) - model_.flux(uL)) / (uR - uL);
    } else {
      linearized_flux = model_.max_eigenvalue(uL);
    }

    if (linearized_flux < 0) {
      return model_.flux(uR);
    } else {
      return model_.flux(uL);
    }
  }

private:
  Model model_;
};

#endif // FVMSCALAR1D_NUMERICAL_FLUX_HPP
