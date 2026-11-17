cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(enkf)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(enkf_libtorch ../enkf_libtorch.cpp)
target_link_libraries(enkf_libtorch "${TORCH_LIBRARIES}")
set_property(TARGET enkf_libtorch PROPERTY CXX_STANDARD 20)

