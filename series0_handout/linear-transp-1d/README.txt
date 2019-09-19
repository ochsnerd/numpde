TEMPLATE CODE FOR SERIES 5
Computational Methods for Engineering Applications II

# Overview

This is the template code for series 5 in the course Computational Methods for Engineering Applications II.

Course homepage:
    http://www2.math.ethz.ch/education/bachelor/lectures/hs2015/other/cmea2

# Requirements
  * A C++ compiler (GCC, Clang and MSVC tested).
  * cmake [optional]
  * EITHER MATLAB or Python (with numpy + matplotlib) for plotting.


# Building with CMake

We strongly recommend building the source using cmake.

## Building with cmake from the command line

In the template_code folder, first make a new build folder:

    mkdir build
    cd build

then issue cmake:

    cmake ..

You can build the whole program using

    make

## Running the executables

If you used cmake, you can now run the executables using

    ./linear_transport
	
## About Eigen

We have bundled the template and solution code with the Eigen library (see http://eigen.tuxfamily.org/). For license information on using Eigen, see the contents of the folder LicenseEigen. *All files under the folder Eigen are subject to the licenses used by Eigen*.
    
