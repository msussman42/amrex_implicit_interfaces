cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(pinn_libtorch)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(pinn_example ../pinn_example.cpp)
target_link_libraries(pinn_example "${TORCH_LIBRARIES}")
set_property(TARGET pinn_example PROPERTY CXX_STANDARD 20)

