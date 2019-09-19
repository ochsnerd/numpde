#pragma once
#include <Eigen/Core>

//! Makes a coordinate transform that maps 
//!
//! e1 to a1
//! e2 to a2
//! 
//! where {e1, e2} is the standard basis for R^2.
template<class Point>
Eigen::Matrix2d makeCoordinateTransform(const Point& a1, const Point& a2) {
    Eigen::Matrix2d coordinateTransform;

    coordinateTransform << a1(0), a2(0),
                           a1(1), a2(1);

    return coordinateTransform;
}
