#pragma once
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

///
/// Converts the argument to string in an OS-independent manner
/// The point here is that NaN and inf may be printed differently, this
/// takes care of that.
///
template <class T>
void printToStream(std::ostream &out, const T &in) {
    if (in != in) { // in is NaN
        out << "NaN";
    } else if (std::isfinite(in)) {
        out << in;
    } else {
        out << "inf";
    }
}

///
/// Writes the contents of the vector 'data' to the textfile filename
/// The output format should be load-able in MATLAB and Numpy using the
/// load-command.
///
template <typename T>
void writeToFile(const std::string &filename, const T &data) {

    std::ofstream file(filename.c_str());
    // You need to check that the stream is good to write to. Even if there is
    // nothing you can do about it, like here.
    assert(file.good());

    // Set highest possible precision, this way we are sure we are
    file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    // Loop over vector and write output to file
    for (size_t i = 0; i < data.size(); ++i) {
        printToStream(file, data[i]);
        file << " ";
    }
    file << std::endl;

    // File closes automatically at end of scope!
}

template <typename T>
void appendMatrixToFile(const std::string &filename, const T &data) {

    std::ofstream file(filename.c_str(), std::fstream::app);
    // You need to check that the stream is good to write to. Even if there is
    // nothing you can do about it, like here.
    assert(file.good());

    // Set highest possible precision, this way we are sure we are
    file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    // Loop over vector and write output to file
    for (int i = 0; i < data.cols(); ++i) {
        for (int j = 0; j < data.rows(); ++j) {
            printToStream(file, data(j, i));
            file << " ";
        }
        file << std::endl;
    }
    file << std::endl;

    // File closes automatically at end of scope!
}
