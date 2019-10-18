# Mandatory information that CMAKE needs
cmake_minimum_required(VERSION 3.12.0)

# Project name
project(ADVNUMCSEFD2D)

# Manage dependencies
include("cmake/GetEigen.cmake")
include("cmake/GetJSON.cmake")

include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

add_executable(2d_transport)
target_compile_features(2d_transport PUBLIC cxx_std_17)
target_sources(2d_transport PRIVATE linear_transport.cpp)
target_link_libraries(2d_transport PRIVATE warning_flags)
target_link_libraries(2d_transport PUBLIC Eigen)
target_link_libraries(2d_transport PUBLIC JSON)


if(CLANG_TIDY_EXE)
  set_target_properties(2d_transport PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()

# So that Visual Studio will define M_PI
target_compile_definitions(2d_transport PUBLIC -D_USE_MATH_DEFINES)
