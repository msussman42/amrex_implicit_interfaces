cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(dynamicMLP)

find_package(Torch REQUIRED)
set(CMAKE_CXX_STANDARD 20)

add_executable(dynamicMLP ../dynamic_MLP.cpp)
target_link_libraries(dynamicMLP "${TORCH_LIBRARIES}")
set_property(TARGET dynamicMLP PROPERTY CXX_STANDARD 20)

