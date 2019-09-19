#include <Eigen/Core>
#include <sstream>
#include "fem_solve.hpp"
#include "writer.hpp"

double f_square(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI *x) * sin(M_PI * y);
}



int main(int argc, char** argv) {
  try {
    Vector u;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi triangles;
    Eigen::MatrixXi tetrahedra;

    if (argc < 2) {
        std::cerr << "Usage: ./fem2d {square,Lshape}_{0..7}" << std::endl;
        return EXIT_FAILURE;
    }
    std::string src = argv[1];

    igl::readMESH(ANCSE_DATA_PATH + src + ".mesh", vertices, tetrahedra, triangles);

    solveFiniteElement(u, vertices, triangles, f_square);

    writeToFile(src + "_values.txt", u);
    writeMatrixToFile(src + "_vertices.txt", vertices);
    writeMatrixToFile(src + "_triangles.txt", triangles);

  }
  catch (std::runtime_error& e) {
    std::cerr << "An error occurred. Error message: " << std::endl;
    std::cerr << "    \"" << e.what() << "\"" << std::endl;
    return EXIT_FAILURE;
  } 
  catch (...) {
    std::cerr << "An unknown error occurred." << std::endl;
    throw ;
  }

  return EXIT_SUCCESS;

}
