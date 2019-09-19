#
# include, download external libraries
#
include(ExternalProject)
#  all external projects should use the same compiler
set(EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX
    "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}" )

find_package(Eigen3)
add_library(Eigen INTERFACE)

if (EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
  target_include_directories(Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
else()
  SET(DOWNLOADING_EIGEN ON)
  #  if not found system wide download
  message("-- Downloading Eigen3")
  ExternalProject_Add(
      EigenProject
      URL http://bitbucket.org/eigen/eigen/get/3.3.7.zip
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
      INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
      CMAKE_ARGS
          ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX}
          -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
  )

  add_dependencies(Eigen EigenProject)
  target_include_directories(Eigen
    SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3
  )

endif()
