# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 3.12.0)


# Project name
PROJECT (2d-poissonlFEM)

set(CMAKE_CXX_STANDARD 11)

find_package (Threads REQUIRED)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # Ignore specific warnings originating in Eigen or IGL
    add_definitions(-Wno-deprecated-register -Wno-return-type-c-linkage)
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  # Ignore specific warnings originating in Eigen or IGL
    add_definitions(-Wno-deprecated-declarations -Wno-ignored-attributes)
endif()

INCLUDE(ExternalProject)

ExternalProject_Add(
    Eigen
    URL http://bitbucket.org/eigen/eigen/get/3.2.7.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(
    IGL
    URL https://github.com/libigl/libigl/archive/master.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IGL
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/Eigen)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/IGL/include)
INCLUDE_DIRECTORIES(".")

ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ADD_DEFINITIONS(-DANCSE_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")

FILE(GLOB SRC *.cpp *.hpp)
ADD_EXECUTABLE(fem2d ${SRC})
TARGET_LINK_LIBRARIES(fem2d ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

ADD_DEPENDENCIES(fem2d IGL Eigen )

ADD_SUBDIRECTORY(unittest)
