cmake_minimum_required (VERSION 3.2.0)



ExternalProject_Add(
    GTest
    URL https://github.com/google/googletest/archive/release-1.8.1.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_source
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
    CMAKE_ARGS "-Dgtest_force_shared_crt=ON"
    INSTALL_COMMAND  ""
    )

# We need the thread package for GTest
FIND_PACKAGE(Threads REQUIRED)



FILE(GLOB_RECURSE SRC *.cpp)

INCLUDE_DIRECTORIES(  ${CMAKE_CURRENT_BINARY_DIR}/gtest_source/googletest/include "..")



LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
  ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary/googlemock/gtest/${CMAKE_BUILD_TYPE}
  ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary/googlemock/gtest)
ADD_EXECUTABLE(unittest ${SRC})
ADD_DEPENDENCIES(unittest GTest IGL Eigen)

TARGET_LINK_LIBRARIES(unittest gtest$<$<AND:$<OR:$<PLATFORM_ID:Windows>,$<AND:$<PLATFORM_ID:Darwin>,$<CXX_COMPILER_ID:AppleClang>>>,$<CONFIG:Debug>>:d> ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
