cmake_minimum_required(VERSION 3.16)

cmake_policy(SET CMP0048 NEW)
project(HDNum-Benchmarks VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(benchmark REQUIRED)

include_directories(../)

macro(package_add_bench BENCHNAME)
    # create an exectuable in which the benchmark will be stored
    add_executable(${BENCHNAME} ${ARGN})
    target_link_libraries(${BENCHNAME} benchmark::benchmark)
endmacro()

package_add_bench(MatrixVectorReference mv.cpp)
package_add_bench(MatrixVectorOperator operator_star.cpp)
package_add_bench(MatrixVectorOperatorDirect ./operator_star_direct.cpp)
