cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(pinn_libtorch)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(wave_4DVAR_noise ../wave_4DVAR_noise.cpp)
target_link_libraries(wave_4DVAR_noise "${TORCH_LIBRARIES}")
set_property(TARGET wave_4DVAR_noise PROPERTY CXX_STANDARD 20)

