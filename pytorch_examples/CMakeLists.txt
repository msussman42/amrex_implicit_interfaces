cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(pinn_libtorch)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(wave_pinnFD ../wave_pinnFD.cpp)
target_link_libraries(wave_pinnFD "${TORCH_LIBRARIES}")
set_property(TARGET wave_pinnFD PROPERTY CXX_STANDARD 20)

