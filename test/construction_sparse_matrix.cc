#include <complex>
#include <initializer_list>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test a non-trivial SparseMatrix class.

using hdnum::SparseMatrix, std::initializer_list;

template <typename T>
class TestSparseMatrix : public ::testing::Test {
public:
    const SparseMatrix<T> sizedConstructed;

    const SparseMatrix<T> fromInitializerQuad;
    const SparseMatrix<T> fromInitializerList;

    using size_type = typename SparseMatrix<T>::size_type;

    constexpr static inline const size_type dimN = 5;
    constexpr static inline const size_type dimM = 5;
    constexpr static inline const size_type value = 5;

    constexpr static inline const initializer_list<initializer_list<T>>
        initializerList = {{0, 1, 2, 3}};
    constexpr static inline const initializer_list<initializer_list<T>>
        initializerListQuad = {
            {0, 1, 2, 3, 4, 5},       {6, 7, 8, 9, 10, 11},
            {12, 13, 14, 15, 16, 17}, {18, 19, 20, 21, 22, 23},
            {24, 25, 26, 27, 28, 29}, {30, 31, 32, 33, 34, 35}};

    TestSparseMatrix()
        : sizedConstructed(dimM, dimN),
          fromInitializerQuad(
              typename SparseMatrix<T>::builder(initializerListQuad).build()),
          fromInitializerList(
              typename SparseMatrix<T>::builder(initializerList).build()) {}
};

using TestTypes = ::testing::Types<int, double, float>;
TYPED_TEST_SUITE(TestSparseMatrix, TestTypes);

TYPED_TEST(TestSparseMatrix, SizeTest) {
    using T = TestSparseMatrix<TypeParam>;
    using size_type = typename T::size_type;

    EXPECT_EQ(T::dimM, this->sizedConstructed.rowsize());
    EXPECT_EQ(T::dimN, this->sizedConstructed.colsize());

    EXPECT_EQ(size_type(6), this->fromInitializerQuad.rowsize());
    EXPECT_EQ(size_type(6), this->fromInitializerQuad.colsize());

    EXPECT_EQ(size_type(1), this->fromInitializerList.rowsize());
    EXPECT_EQ(size_type(4), this->fromInitializerList.colsize());
}

TYPED_TEST(TestSparseMatrix, ValueIndexTest) {
    using T = TestSparseMatrix<TypeParam>;
    using size_type = typename T::size_type;

    for (auto i = size_type(0); i < this->sizedConstructed.rowsize(); i++)
        for (auto j = size_type(0); j < this->sizedConstructed.colsize(); j++)
            EXPECT_EQ(TypeParam(0), this->sizedConstructed(i, j));

    for (auto i = size_type(0); i < this->fromInitializerList.rowsize(); i++)
        for (auto j = size_type(0); j < this->fromInitializerList.colsize();
             j++)
            EXPECT_EQ(TypeParam(i * this->fromInitializerList.rowsize() + j),
                      this->fromInitializerList(i, j));

    for (auto i = size_type(0); i < this->fromInitializerQuad.rowsize(); i++)
        for (auto j = size_type(0); j < this->fromInitializerQuad.colsize();
             j++)
            EXPECT_EQ(TypeParam(i * this->fromInitializerQuad.rowsize() + j),
                      this->fromInitializerQuad(i, j));
}

TYPED_TEST(TestSparseMatrix, ConversionTest) {
    using size_type = typename SparseMatrix<TypeParam>::size_type;
    auto denseMatrix = hdnum::DenseMatrix(this->fromInitializerQuad);

    EXPECT_EQ(denseMatrix.rowsize(), this->fromInitializerQuad.rowsize());
    EXPECT_EQ(denseMatrix.colsize(), this->fromInitializerQuad.colsize());

    for (size_type i = 0; i < denseMatrix.rowsize(); i++) {
        for (size_type j = 0; j < denseMatrix.colsize(); j++) {
            EXPECT_EQ(denseMatrix[i][j], this->fromInitializerQuad(i, j));
        }
    }
}

}  // namespace
