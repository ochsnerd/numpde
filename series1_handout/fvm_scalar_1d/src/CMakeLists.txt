add_subdirectory(ancse)

target_sources(fvm_scalar
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fvm_scalar.cpp
)
