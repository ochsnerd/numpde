#pragma once
#include <Eigen/Core>
#include <functional>


//! Approximates the integral
//! 
//! \f[\int_K f(x, y) \; dV\f]
//!
//! where K is the triangle spanned by (0,0), (1,0) and (0, 1).
inline double integrate(const std::function<double(double, double)>& f)
{
    Eigen::VectorXd weights(7);
    weights << 9. / 80.,
        (155 + sqrt(15)) / 2400,
        (155 + sqrt(15)) / 2400,
        (155 + sqrt(15)) / 2400,
        (155 - sqrt(15)) / 2400,
        (155 - sqrt(15)) / 2400,
        (155 - sqrt(15)) / 2400;


    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> points(7, 2);
    points << 1.0 / 3., 1. / 3., 
        (6 + sqrt(15)) / 21.,  (6 + sqrt(15)) / 21.,
        (9 - 2 * sqrt(15)) / 21,   (6 + sqrt(15)) / 21,
        (6 + sqrt(15)) / 21,  (9 - 2 * sqrt(15)) / 21,
        (6 - sqrt(15)) / 21, (9 + 2 * sqrt(15)) / 21,
        (9 + 2 * sqrt(15)) / 21,   (6 - sqrt(15)) / 21,
        (6 - sqrt(15)) / 21,  (6 - sqrt(15)) / 21;
    
    double integral = 0.0;
    for (int i = 0; i < points.rows(); ++i) {
        integral += weights(i) * f(points(i, 0), points(i, 1));
    }

    return integral;
}
