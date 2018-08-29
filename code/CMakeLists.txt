cmake_minimum_required(VERSION 3.5)

project(proyecto1 VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory (src)

enable_testing ()
add_subdirectory (test)
add_subdirectory(benchmarks)
