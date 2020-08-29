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
    const SparseMatrix<T> sizedConstructedZero;
    const SparseMatrix<T> sizedConstructedValue;

    const SparseMatrix<T> fromInitializerListQuad;
    const SparseMatrix<T> fromInitializerList;

    const SparseMatrix<T> fromDenseMatrix;

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

    static inline const auto someDenseMatrix =
        hdnum::DenseMatrix<T>(initializerListQuad);

    TestSparseMatrix()
        : sizedConstructed(dimM, dimN), sizedConstructedZero(dimM, dimN, T(0)),
          sizedConstructedValue(dimM, dimN, T(value)),
          fromInitializerListQuad(initializerListQuad),
          fromInitializerList(initializerList),
          fromDenseMatrix(someDenseMatrix) {}
};

using TestTypes = ::testing::Types<int, double, float, std::complex<int>,
                                   std::complex<double>, std::complex<float>>;
TYPED_TEST_SUITE(TestSparseMatrix, TestTypes);

TYPED_TEST(TestSparseMatrix, SizeTest) {
    using T = TestSparseMatrix<TypeParam>;
    using size_type = typename T::size_type;

    ASSERT_EQ(T::dimM, this->sizedConstructed.rowsize());
    ASSERT_EQ(T::dimN, this->sizedConstructed.colsize());

    ASSERT_EQ(T::dimM, this->sizedConstructedZero.rowsize());
    ASSERT_EQ(T::dimN, this->sizedConstructedZero.colsize());

    ASSERT_EQ(T::dimM, this->sizedConstructedValue.rowsize());
    ASSERT_EQ(T::dimN, this->sizedConstructedValue.colsize());

    ASSERT_EQ(size_type(4), this->fromInitializerListQuad.rowsize());
    ASSERT_EQ(size_type(4), this->fromInitializerListQuad.colsize());

    ASSERT_EQ(size_type(4), this->fromInitializerList.rowsize());
    ASSERT_EQ(size_type(1), this->fromInitializerList.colsize());

    ASSERT_EQ(T::someDenseMatrix.rowsize(), this->fromDenseMatrix.rowsize());
    ASSERT_EQ(T::someDenseMatrix.colsize(), this->fromDenseMatrix.colsize());
}

TYPED_TEST(TestSparseMatrix, ValueIndexTest) {
    using T = TestSparseMatrix<TypeParam>;
    using size_type = typename T::size_type;

    for (auto i = size_type(0); this->sizedConstructed.rowsize(); i++)
        for (auto j = size_type(0); this->sizedConstructed.colsize(); j++)
            ASSERT_EQ(TypeParam(0), this->sizedConstructed(i, j));

    for (auto i = size_type(0); this->sizedConstructedZero.rowsize(); i++)
        for (auto j = size_type(0); this->sizedConstructedZero.colsize(); j++)
            ASSERT_EQ(TypeParam(0), this->sizedConstructedZero(i, j));

    for (auto i = size_type(0); this->sizedConstructedValue.rowsize(); i++)
        for (auto j = size_type(0); this->sizedConstructedValue.colsize(); j++)
            ASSERT_EQ(TypeParam(T::value), this->sizedConstructedValue(i, j));

    for (auto i = size_type(0); this->fromInitializerList.rowsize(); i++)
        for (auto j = size_type(0); this->fromInitializerList.colsize(); j++)
            ASSERT_EQ(TypeParam(i * this->fromInitializerList.rowsize() + j),
                      this->fromInitializerList(i, j));

    for (auto i = size_type(0); this->fromInitializerListQuad.rowsize(); i++)
        for (auto j = size_type(0); this->fromInitializerListQuad.colsize();
             j++)
            ASSERT_EQ(
                TypeParam(i * this->fromInitializerListQuad.rowsize() + j),
                this->fromInitializerListQuad(i, j));

    for (auto i = size_type(0); this->fromDenseMatrix.rowsize(); i++)
        for (auto j = size_type(0); this->fromDenseMatrix.colsize(); j++)
            ASSERT_EQ(T::someDenseMatrix(i, j), this->fromDenseMatrix(i, j));
}

/* TYPED_TEST(TestSparseMatrix, ValueConstIteratorTest) { */
/*     for (const auto& row : this->sizedConstructed) */
/*         for (const auto& element : row) ASSERT_EQ(TypeParam(0), element); */

/*     for (const auto& row : this->sizedConstructedZero) */
/*         for (const auto& element : row) ASSERT_EQ(TypeParam(0), element); */

/*     for (const auto& row : this->sizedConstructedValue) */
/*         for (const auto& element : row) */
/*             ASSERT_EQ(TypeParam(this->value), element); */
/* } */

/* TYPED_TEST(TestSparseMatrix, ValueIteratorTest) { */
/*     for (auto& row : this->sizedConstructed) */
/*         for (auto& element : row) ASSERT_EQ(TypeParam(0), element); */

/*     for (auto& row : this->sizedConstructedZero) */
/*         for (auto& element : row) ASSERT_EQ(TypeParam(0), element); */

/*     for (auto& row : this->sizedConstructedValue) */
/*         for (auto& element : row) ASSERT_EQ(TypeParam(this->value), element);
 */
/* } */

}  // namespace
