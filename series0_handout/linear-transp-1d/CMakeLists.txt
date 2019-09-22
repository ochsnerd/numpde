# Mandatory information that CMAKE needs
cmake_minimum_required(VERSION 3.12.0)

# Project name
project(ADVNUMCSEFD1D)

# Manage dependencies
# include("cmake/GetEigen.cmake")

include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

add_executable(linear_transport)
set_target_properties(linear_transport PROPERTIES CXX_STANDARD 17)
target_sources(linear_transport PRIVATE linear_transport.cpp)
target_link_libraries(linear_transport PRIVATE warning_flags)
target_include_directories(linear_transport PUBLIC $(CMAKE_SOURCE_DIR))
# target_link_libraries(linear_transport PUBLIC Eigen)

# So that Visual Studio will define M_PI
target_compile_definitions(linear_transport PUBLIC -D_USE_MATH_DEFINES)
