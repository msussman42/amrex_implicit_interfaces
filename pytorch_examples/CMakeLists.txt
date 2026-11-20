cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(enkf)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(advect_4DVAR ../advect_4DVAR.cpp)
target_link_libraries(advect_4DVAR "${TORCH_LIBRARIES}")
set_property(TARGET advect_4DVAR PROPERTY CXX_STANDARD 20)

