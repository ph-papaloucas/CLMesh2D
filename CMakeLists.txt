cmake_minimum_required(VERSION 3.16)
project(PDESolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
include_directories(include)

# Collect all source files
file(GLOB_RECURSE SOURCES src/*.cc)

# Find all example files
file(GLOB EXAMPLES examples/*.cc)
file (GLOB TESTS tests/*.cc)

# Build an executable for each example
foreach(example ${EXAMPLES})
    get_filename_component(example_name ${example} NAME_WE)
    add_executable(${example_name} ${example} ${SOURCES})
endforeach()

foreach (test ${TESTS})
    get_filename_component(test_name ${test} NAME_WE)
    add_executable(${test_name} ${test} ${SOURCES})
endforeach()
