/* #include <complex> */

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test a non-trivial SparseMatrix class.

using hdnum::SparseMatrix, std::initializer_list;

template <typename T>
class TestSparseMatrixBuilder : public ::testing::Test {};

/* using TestTypes = ::testing::Types<int, double, float, std::complex<int>, */
/*                                    std::complex<double>,
 * std::complex<float>>; */

using TestTypes = ::testing::Types<int, double>;
TYPED_TEST_SUITE(TestSparseMatrixBuilder, TestTypes);

TYPED_TEST(TestSparseMatrixBuilder, SizeTest) {
    /* using T = TestSparseMatrixBuilder<TypeParam>; */
    auto builder = typename SparseMatrix<TypeParam>::builder();
}

}  // namespace
