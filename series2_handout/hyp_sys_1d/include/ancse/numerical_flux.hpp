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

class RoeFlux : public NumericalFlux<RoeFlux> {
public:
  using Vector = typename NumericalFlux<RoeFlux>::Vector;

  explicit RoeFlux(const std::shared_ptr<Model>& model) : model_{model} {}

  virtual Vector compute_flux(const Vector& uL, const Vector& uR) const override {
    auto flux_average = 0.5 * (model_->flux(uR) + model_->flux(uL));
    auto diffusion = 0.5 * eigenbasis_(uL, uR) * abs_(eigenvalues_(uL, uR)) * inv_eigenbasis_(uL, uR) * (uR - uL);
    return  flux_average - diffusion;
  }

private:
  using Matrix = Eigen::MatrixXd;

  Matrix eigenbasis_(const Vector& uL, const Vector& uR) const {
    double v_bar = roe_average(uL, uR, v(uL), v(uR));
    double h_bar = roe_average(uL, uR, H(uL), H(uR));
    double c_bar = std::sqrt((1 - gamma) * (h_bar - 0.5 * v_bar * v_bar));

    Matrix R(3,3);
    R <<
      1,                     1,                  1,
      v_bar - c_bar,         v_bar,              v_bar + c_bar,
      h_bar - v_bar * c_bar, .5 * v_bar * v_bar, h_bar + v_bar * c_bar;
    return R;
  }

  Matrix eigenvalues_(const Vector& uL, const Vector& uR) const {
    double v_bar = roe_average(uL, uR, v(uL), v(uR));
    double h_bar = roe_average(uL, uR, H(uL), H(uR));
    double c_bar = std::sqrt((1 - gamma) * (h_bar - 0.5 * v_bar * v_bar));
    return Eigen::Vector3d(v_bar - c_bar, v_bar, v_bar + c_bar).asDiagonal();
  }

  Matrix inv_eigenbasis_(const Vector& uL, const Vector& uR) const {
    // plz invert eigenbasis_(uL, uR)
    // make sure to handle normalization correctly
    return Matrix{};
  }

  Matrix abs_(const Matrix& lambda) const {
    // Could implement Hartens entropy fix here for example
    return abs(lambda.array());
  }

  double roe_average(const Vector& uL, const Vector& uR, double xL, double xR) const {
    double sqrt_rhoL = std::sqrt(uL[0]);
    double sqrt_rhoR = std::sqrt(uR[0]);

    return (sqrt_rhoL * xL + sqrt_rhoR * xR) / (sqrt_rhoL + sqrt_rhoR);
  }

  // Dilemma: This should semantically be handled by the model. However we work
  // with a pointer to the base class, so we'd have to implement these functions
  // also for other models (Burgers) where they dont make sense.
  // A solution would be to have a RoeModel thats an adaptor to model and implements
  // methods to get the eigendecomposition of the Roe matrix (or a Roe-Matrix class).
  double rho(const Vector& u) const {
    return u[0];
  }

  double v(const Vector& u) const {
    return u[1] / u[0];
  }

  double H(const Vector& u) const {
    double E = u[2];
    double p = (E - .5 * rho(u) * v(u) * v(u)) * (gamma - 1);
    return (E + p) / rho(u);
  }

  double gamma = 5./3;  // this is horrible
  std::shared_ptr<Model> model_;
};
                                       
// template<class Diffusion>
// class ApproximateFlux : public NumericalFlux<ApproximateFlux>
// // Fluxes of the form "flux-average + diffusion"
// public:
//   virtual ~ApproximateFlux() {};

// does this make sense? almost surely overengineered

#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
