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

TYPED_TEST(TestSparseMatrixBuilder, EqualityComparison) {
    // Example from: https://de.wikipedia.org/wiki/Compressed_Row_Storage
    auto builder1 = typename SparseMatrix<TypeParam>::builder(4, 5);
    auto builder2 = typename SparseMatrix<TypeParam>::builder(4, 5);

    EXPECT_EQ(builder1.rowsize(), 4);
    EXPECT_EQ(builder1.colsize(), 5);

    EXPECT_EQ(builder2.rowsize(), 4);
    EXPECT_EQ(builder2.colsize(), 5);

    builder1.addEntry(3, 2, 11);
    builder1.addEntry(0, 0, 10);
    builder1.addEntry(1, 2, 11);
    builder1.addEntry(0, 3, 12);
    builder1.addEntry(2, 1, 16);
    builder1.addEntry(1, 4, 13);
    builder1.addEntry(3, 4, 13);

    builder2.addEntry(3, 2, 11);
    builder2.addEntry(0, 0, 10);
    builder2.addEntry(1, 2, 11);
    builder2.addEntry(0, 3, 12);
    builder2.addEntry(2, 1, 16);
    builder2.addEntry(1, 4, 13);
    builder2.addEntry(3, 4, 13);

    EXPECT_EQ(builder1, builder2);

    auto M1 = builder1.build();
    auto M2 = builder2.build();

    EXPECT_EQ(M1, M2);

    EXPECT_EQ(M1.rowsize(), 4);
    EXPECT_EQ(M1.colsize(), 5);
    EXPECT_EQ(M2.rowsize(), 4);
    EXPECT_EQ(M2.colsize(), 5);

    builder2.addEntry(3, 3, 13);
    M2 = builder2.build();

    EXPECT_NE(builder1, builder2);
    EXPECT_NE(M1, M2);
}

}  // namespace
