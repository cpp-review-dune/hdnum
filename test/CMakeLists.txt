cmake_minimum_required(VERSION 3.16)

cmake_policy(SET CMP0048 NEW)
project(HDNum-Tests VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(PkgConfig)
pkg_search_module(GTEST REQUIRED gtest_main)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")

include_directories(../)

macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${TESTNAME} ${GTEST_LDFLAGS})

    target_compile_options(${TESTNAME} PUBLIC ${GTEST_CFLAGS})
    add_test("${TESTNAME}Test" ${TESTNAME})
endmacro()

include(CTest)
package_add_test(SparseMatrixDefault default_sparse_matrix.cc)
package_add_test(SparseMatrixSpecialMatricies special_matricies_sparse_matrix.cc)
package_add_test(SparseMatrixConstruction construction_sparse_matrix.cc)
package_add_test(SparseMatrixOperators operators_sparse_matrix.cc)