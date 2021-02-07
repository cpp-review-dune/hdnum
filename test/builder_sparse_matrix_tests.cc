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

TYPED_TEST(TestSparseMatrixBuilder, DefaultConstructor) {
    auto builder = typename SparseMatrix<TypeParam>::builder();
    EXPECT_EQ(builder.rowsize(), 0);
    EXPECT_EQ(builder.colsize(), 0);
}

TYPED_TEST(TestSparseMatrixBuilder, SizedConstructor) {
    auto builder = typename SparseMatrix<TypeParam>::builder(4, 5);
    EXPECT_EQ(builder.rowsize(), 4);
    EXPECT_EQ(builder.colsize(), 5);

    builder.setNumRows(20);
    builder.setNumCols(42);
    EXPECT_EQ(builder.rowsize(), 20);
    EXPECT_EQ(builder.colsize(), 42);
}

TYPED_TEST(TestSparseMatrixBuilder, AddExistingElements) {
    auto builder = typename SparseMatrix<TypeParam>::builder(4, 5);
    EXPECT_EQ(builder.rowsize(), 4);
    EXPECT_EQ(builder.colsize(), 5);

    builder.addEntry(1, 0, 20);
    builder.addEntry(1, 0, 20);
    builder.addEntry(1, 0, 20);

    builder.addEntry(1, 0);
    builder.addEntry(1, 0);
    builder.addEntry(1, 0);
}

TYPED_TEST(TestSparseMatrixBuilder, AddElements) {
    auto builder = typename SparseMatrix<TypeParam>::builder(4, 5);
    EXPECT_EQ(builder.rowsize(), 4);
    EXPECT_EQ(builder.colsize(), 5);

    builder.addEntry(0, 1, 1);
    builder.addEntry(0, 0, 0);
    builder.addEntry(1, 1, 2);
}

TYPED_TEST(TestSparseMatrixBuilder, WikipediaTestCase) {
    // Example from: https://de.wikipedia.org/wiki/Compressed_Row_Storage
    auto builder = typename SparseMatrix<TypeParam>::builder(4, 5);
    EXPECT_EQ(builder.rowsize(), 4);
    EXPECT_EQ(builder.colsize(), 5);

    builder.addEntry(3, 2, 11);
    builder.addEntry(0, 0, 10);
    builder.addEntry(1, 2, 11);
    builder.addEntry(0, 3, 12);
    builder.addEntry(2, 1, 16);
    builder.addEntry(1, 4, 13);
    builder.addEntry(3, 4, 13);

    auto A = builder.build();
    EXPECT_EQ(A.rowsize(), 4);
    EXPECT_EQ(A.colsize(), 5);

    // TODO: Add checks for the contents of A!
}

}  // namespace
