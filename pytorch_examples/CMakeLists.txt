cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(pinn_libtorch)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(wave_pinn ../wave_pinn.cpp)
target_link_libraries(wave_pinn "${TORCH_LIBRARIES}")
set_property(TARGET wave_pinn PROPERTY CXX_STANDARD 20)

