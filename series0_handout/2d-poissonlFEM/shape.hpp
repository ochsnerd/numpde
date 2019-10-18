#pragma once
#include<stdexcept>

//! The shape function (on the reference element)
//! 
//! We have three shape functions.
//! 
//! lambda(0, x, y) should be 1 in the point (0,0) and zero in (1,0) and (0,1)
//! lambda(1, x, y) should be 1 in the point (1,0) and zero in (0,0) and (0,1)
//! lambda(2, x, y) should be 1 in the point (0,1) and zero in (0,0) and (1,0)
//!
//! @param i integer between 0 and 2 (inclusive). Decides which shape function to return.
//! @param x x coordinate in the reference element.
//! @param y y coordinate in the reference element.
inline double lambda(int i, double x, double y) {
  switch(i) {
  case 0:
    return 1 - x - y;
  case 1:
    return x;
  case 2:
    return y;
  default:
    throw std::invalid_argument("i has to be in {0,1,2}"); break;
  }
}
