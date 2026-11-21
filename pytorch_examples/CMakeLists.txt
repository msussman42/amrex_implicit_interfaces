cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(enkf)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(advect_ENKF ../advect_ENKF.cpp)
target_link_libraries(advect_ENKF "${TORCH_LIBRARIES}")
set_property(TARGET advect_ENKF PROPERTY CXX_STANDARD 20)

