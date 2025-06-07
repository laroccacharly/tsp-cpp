cmake_minimum_required(VERSION 3.20)
project(tsp-cpp)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(fmt CONFIG REQUIRED)
find_package(highs CONFIG REQUIRED)

add_executable(tsp-cpp main.cpp)
target_link_libraries(tsp-cpp PRIVATE fmt::fmt highs::highs)