cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(data_assimilation)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(advect_compare_ENKF_kalman ../advect_compare_ENKF_kalman.cpp)
target_link_libraries(advect_compare_ENKF_kalman "${TORCH_LIBRARIES}")
set_property(TARGET advect_compare_ENKF_kalman PROPERTY CXX_STANDARD 20)

