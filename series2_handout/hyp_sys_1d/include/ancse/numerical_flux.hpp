#ifndef HYPSYS1D_NUMERICAL_FLUX_HPP
#define HYPSYS1D_NUMERICAL_FLUX_HPP

#include <memory>
#include <iostream>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model.
 * It is also unconditionally a bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access
    //       to the following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const std::shared_ptr<Model> &model)
        : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    std::shared_ptr<Model> model;
};

template<class ConcreteFlux>
class NumericalFlux {
  // Interface for numerical fluxes

  // Using CRTP to get compile time polymorphism
public:
  using Vector = Eigen::VectorXd;

  virtual ~NumericalFlux() {};

  virtual Vector compute_flux(const Vector&, const Vector&) const = 0;

  Vector operator() (const Vector& uL, const Vector& uR) const {
    return static_cast<ConcreteFlux const*>(this)->compute_flux(uL, uR);
  }

  Vector operator() (const std::pair<Vector,Vector>& u) const {
    return static_cast<ConcreteFlux const*>(this)->compute_flux(u.first, u.second);
  }
};


class RusanovFlux : public NumericalFlux<RusanovFlux> {
  // Mishra Hyperbolic PDES 4.2.4
public:
  using Vector = typename NumericalFlux<RusanovFlux>::Vector;

  explicit RusanovFlux(const std::shared_ptr<Model>& model) : model_{model} {}

  Vector compute_flux(const Vector& uL, const Vector& uR) const override
  {
    double speed = std::max(std::abs(model_->max_eigenvalue(uL)),
                            std::abs(model_->max_eigenvalue(uR)));

    return 0.5 * (model_->flux(uR) + model_->flux(uL)) - 0.5 * speed * (uR - uL);
  }

private:
  std::shared_ptr<Model> model_;
};

class LaxFriedrichsFlux : public NumericalFlux<LaxFriedrichsFlux> {
  // LN 4.2.3
public:
  using Vector = typename NumericalFlux<LaxFriedrichsFlux>::Vector;

  LaxFriedrichsFlux(const std::shared_ptr<Model>& model,
                    const Grid& grid,
                    const std::shared_ptr<SimulationTime>& time) :
    model_{model},
    speed_{grid.dx / time->dt}
  {}

  virtual Vector compute_flux(const Vector& uL, const Vector& uR) const override {
    return 0.5 * (model_->flux(uR) + model_->flux(uL)) - 0.5 * speed_ * (uR -uL);
  }

private:
  std::shared_ptr<Model> model_;
  double speed_;
};
// template<class Diffusion>
// class ApproximateFlux : public NumericalFlux<ApproximateFlux>
// // Fluxes of the form "flux-average + diffusion"
// public:
//   virtual ~ApproximateFlux() {};

// does this make sense? almost surely overengineered

#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
