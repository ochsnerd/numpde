# Mandatory information that CMAKE needs
cmake_minimum_required(VERSION 3.12.0)

# Project name
project(FVMScalar1D)

# Manage dependencies
include("cmake/GetEigen.cmake")
include("cmake/GetFMT.cmake")
include("cmake/GetJSON.cmake")

include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

add_library(fvm_scalar_obj OBJECT)
target_compile_features(fvm_scalar_obj PUBLIC cxx_std_17)
target_include_directories(fvm_scalar_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(fvm_scalar_obj PUBLIC Eigen)
target_link_libraries(fvm_scalar_obj PUBLIC FMT)
target_link_libraries(fvm_scalar_obj PUBLIC JSON)
target_link_libraries(fvm_scalar_obj PRIVATE warning_flags)
target_compile_definitions(fvm_scalar_obj PUBLIC -D_USE_MATH_DEFINES)

if(CLANG_TIDY_EXE)
  set_target_properties(fvm_scalar_obj PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()

add_executable(fvm_scalar)
target_link_libraries(fvm_scalar PRIVATE fvm_scalar_obj)

# Testing
include(cmake/GetGTest.cmake)

add_executable(unit_tests)
target_link_libraries(unit_tests fvm_scalar_obj)
target_link_libraries(unit_tests GTest)

add_subdirectory(src)
add_subdirectory(tests)
