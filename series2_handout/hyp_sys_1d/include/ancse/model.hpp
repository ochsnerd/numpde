#ifndef HYPSYS1D_MODEL_HPP
#define HYPSYS1D_MODEL_HPP

#include <cmath>
#include <memory>

#include <Eigen/Dense>
#include <ancse/config.hpp>

/// Interface for implementing different models,
/// eg. Euler equations, Shallow-water equations
///
/// Add more functions to this interface if needed.
class Model {
  public:
  using Vector = Eigen::VectorXd;
  using Matrix = Eigen::MatrixXd;
  virtual ~Model() = default;

  virtual Vector flux(const Vector &u) const = 0;
  virtual Vector eigenvalues(const Vector &u) const = 0;
  virtual Matrix eigenvectors(const Vector &u) const = 0;
  virtual double max_eigenvalue(const Vector &u) const = 0;

  virtual Vector cons_to_prim(const Vector &u) const = 0;
  virtual Vector prim_to_cons(const Vector &u) const = 0;
  virtual void cons_to_prim(const Matrix& u_cons, Matrix& u_prim) const = 0;
  virtual void prim_to_cons(const Matrix& u_prim, Matrix& u_cons) const = 0;

  virtual int get_nvars() const = 0;
  virtual std::string get_name() const = 0;
};

class Burgers : public Model {
  public:

    Eigen::VectorXd flux(const Eigen::VectorXd &u) const override
    {
        Eigen::VectorXd f(n_vars);
        f(0) = 0.5*u(0)*u(0);

        return f;
    }

    Eigen::VectorXd eigenvalues(const Eigen::VectorXd &u) const override
    {
        Eigen::VectorXd eigvals(n_vars);
        eigvals(0) = u(0);

        return eigvals;
    }

    Eigen::MatrixXd eigenvectors(const Eigen::VectorXd &) const override
    {
        Eigen::MatrixXd eigvecs(n_vars, n_vars);
        eigvecs (0,0) = 1;

        return eigvecs;
    }

    double max_eigenvalue(const Eigen::VectorXd &u) const override {
        return (eigenvalues(u).cwiseAbs()).maxCoeff();
    }

    Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u) const override {
        return u;
    }

    void cons_to_prim(const Eigen::MatrixXd& u_cons, Eigen::MatrixXd& u_prim) const override {
      u_prim = u_cons;
    }

    Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u) const override {
        return u;
    }

    void prim_to_cons(const Eigen::MatrixXd& u_prim, Eigen::MatrixXd& u_cons) const override {
      u_cons = u_prim;
    }

    int get_nvars() const override
    {
        return n_vars;
    }

    std::string get_name() const override
    {
        return "burgers";
    }

  private:
    int n_vars = 1;
};

/// Euler equations
class Euler : public Model {
  // u(0) rho          density
  // u(1) v            velocity
  // u(2) p            pressure
public:
  using Vector = Model::Vector;
  using Matrix = Model::Matrix;

  Vector flux(const Vector &u) const override {
    Vector f(3);
    f <<
      rho(u) * v(u),
      rho(u) * v(u) * v(u) + p(u),
      (E(u) + p(u)) * v(u);
    return f;
  }

  Vector eigenvalues(const Vector &u) const override {
    Vector ev(3);
    ev <<
      v(u) - c(u),
      v(u),
      v(u) + c(u);
    return ev;
  }

  Matrix eigenvectors(const Vector &u) const override {
    Matrix R(3,3);
    R.col(0) << 1.0, v(u) - c(u), H(u) - v(u)*c(u);
    R.col(1) << 1.0, v(u), 0.5 * v(u) * v(u);
    R.col(2) << 1.0, v(u) + c(u), H(u) + v(u)*c(u);
    return R;
  }

  double max_eigenvalue(const Vector &u) const override {
    return eigenvalues(u)(2);
  }

  // Primitive variables:
  // rho, v, p
  Vector cons_to_prim(const Vector &u_cons) const override {
    assert(u_cons.size() == n_vars);
    assert(u_cons(0) > 0);
    assert(u_cons(2) >= 0);

    Vector u_prim(3);
    // u_prim <<
    //   u_cons(0),
    //   u_cons(1) / u_cons(0),
    //   (u_cons(2) - 0.5 * u_cons(1) * u_cons(1) / u_cons(0)) * (gamma - 1);
    u_prim <<
      rho(u_cons),
      v(u_cons),
      p(u_cons);

    return u_prim;
  }

  void cons_to_prim(const Matrix& u_cons, Matrix& u_prim) const override {
    u_prim.resizeLike(u_cons);

    for (int i = 0; i < u_cons.cols(); ++i) {
      u_prim.col(i) = cons_to_prim(u_cons.col(i));
    }
  }

  // Conservative variables:
  // rho, m, E
  Vector prim_to_cons(const Vector &u_prim) const override {
    assert(u_prim.size() == n_vars);
    assert(u_prim(0) > 0);
    assert(u_prim(2) >= 0);

    Vector u_cons(3);
    // u_cons <<
    //   rho(u_prim),
    //   m(u_prim),
    //   E(u_prim);
    u_cons <<
      u_prim(0),
      u_prim(1) * u_prim(0),
      u_prim(2) / (gamma - 1) + .5 * u_prim(1) * u_prim(1) * u_prim(0);

    return u_cons;
  }

  void prim_to_cons(const Matrix& u_prim, Matrix& u_cons) const override {
    u_cons.resizeLike(u_prim);

    for (int i = 0; i < u_prim.cols(); ++i) {
      u_cons.col(i) = prim_to_cons(u_prim.col(i));
    }
  }

  void set_gamma(double gamma_)
  {
    gamma = gamma_;
  }

  double get_gamma() const
  {
    return gamma;
  }

  int get_nvars() const override
  {
    return n_vars;
  }

  std::string get_name() const override
  {
    return "euler";
  }

private:
  // making formulas more readable
  // density
  double rho(const Vector& u) const {
    return u(0);
  }
  // momentum
  double m(const Vector& u) const {
    return u(1);
  }
  // Energy
  double E(const Vector& u) const {
    return u(2);
  }
  // velocity
  double v(const Vector& u) const {
    return m(u) / rho(u);
  }
  // pressure
  double p(const Vector& u) const {
    return (E(u) - .5 * m(u) * m(u) / rho(u)) * (gamma - 1);
  }
  // soundspeed
  double c(const Vector& u) const {
    return std::sqrt(gamma * p(u) / rho(u));
  }
  // Enthalpy
  double H(const Vector u) const {
    return (E(u) + p(u)) / rho(u);  // also only fine when strictly hyperbolic
  }

  int n_vars = 3;
  double gamma = 5./3.;
};


std::shared_ptr<Model> make_model (const nlohmann::json &config);

#endif // HYPSYS1D_MODEL_HPP
