cmake_minimum_required(VERSION 3.16)  # CMake version check

project(lanczos VERSION 1.0) 
set(CMAKE_CXX_STANDARD 17)            # Enable c++17 standard
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # Export compile commands

#include_directories(../../)

# Add main.cpp file of the project root directory as a source file
set(SOURCE_FILES lanczos-example.cc)


# Add executable target with source files listed in SOURCE_FILES variable
add_executable(simple_example ${SOURCE_FILES})