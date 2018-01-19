cmake_minimum_required(VERSION 3.9)
project(asimon C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

add_executable(asimon main.c vector.c vector.h)

find_package(GSL REQUIRED)
target_link_libraries(asimon GSL::gsl GSL::gslcblas)
