#include <complex>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test a non-trivial SparseMatrix class.

using hdnum::SparseMatrix;

template <typename T>
class TestSparseMatrix : public ::testing::Test {
public:
    const SparseMatrix<T> sizedConstructed;
    const SparseMatrix<T> sizedConstructedZero;
    const SparseMatrix<T> sizedConstructedValue;

    const SparseMatrix<T> initializerListQuad;
    const SparseMatrix<T> initializerList;

    using size_type = typename SparseMatrix<T>::size_type;

    constexpr static inline const size_type dimN = 5;
    constexpr static inline const size_type dimM = 5;
    constexpr static inline const size_type value = 5;

    TestSparseMatrix()
        : sizedConstructed(dimM, dimN), sizedConstructedZero(dimM, dimN, T(0)),
          sizedConstructedValue(dimM, dimN, T(value)),
          initializerListQuad({{0, 1, 2, 3, 4, 5},
                               {6, 7, 8, 9, 10, 11},
                               {12, 13, 14, 15, 16, 17},
                               {18, 19, 20, 21, 22, 23},
                               {24, 25, 26, 27, 28, 29},
                               {30, 31, 32, 33, 34, 35}}),
          initializerList({{0, 1, 2, 3}}) {}
};

using TestTypes = ::testing::Types<int, double, float, std::complex<int>,
                                   std::complex<double>, std::complex<float>>;
TYPED_TEST_SUITE(TestSparseMatrix, TestTypes);

TYPED_TEST(TestSparseMatrix, SizeTest) {
    using T = TestSparseMatrix<TypeParam>;
    using size_type = typename T::size_type;

    EXPECT_EQ(T::dimM, this->sizedConstructed.rowsize());
    EXPECT_EQ(T::dimN, this->sizedConstructed.colsize());

    EXPECT_EQ(T::dimM, this->sizedConstructedZero.rowsize());
    EXPECT_EQ(T::dimN, this->sizedConstructedZero.colsize());

    EXPECT_EQ(T::dimM, this->sizedConstructedValue.rowsize());
    EXPECT_EQ(T::dimN, this->sizedConstructedValue.colsize());

    EXPECT_EQ(size_type(4), this->initializerListQuad.rowsize());
    EXPECT_EQ(size_type(4), this->initializerListQuad.colsize());

    EXPECT_EQ(size_type(4), this->initializerList.rowsize());
    EXPECT_EQ(size_type(1), this->initializerList.colsize());
}

TYPED_TEST(TestSparseMatrix, ValueIndexTest) {
    using T = TestSparseMatrix<TypeParam>;
    using size_type = typename T::size_type;

    for (auto i = size_type(0); this->sizedConstructed.rowsize(); i++)
        for (auto j = size_type(0); this->sizedConstructed.colsize(); j++)
            EXPECT_EQ(TypeParam(0), this->sizedConstructed(i, j));

    for (auto i = size_type(0); this->sizedConstructedZero.rowsize(); i++)
        for (auto j = size_type(0); this->sizedConstructedZero.colsize(); j++)
            EXPECT_EQ(TypeParam(0), this->sizedConstructedZero(i, j));

    for (auto i = size_type(0); this->sizedConstructedValue.rowsize(); i++)
        for (auto j = size_type(0); this->sizedConstructedValue.colsize(); j++)
            EXPECT_EQ(TypeParam(T::value), this->sizedConstructedValue(i, j));

    for (auto i = size_type(0); this->initializerList.rowsize(); i++)
        for (auto j = size_type(0); this->initializerList.colsize(); j++)
            EXPECT_EQ(TypeParam(i * this->initializerList.rowsize() + j),
                      this->initializerList(i, j));

    for (auto i = size_type(0); this->initializerListQuad.rowsize(); i++)
        for (auto j = size_type(0); this->initializerListQuad.colsize(); j++)
            EXPECT_EQ(TypeParam(i * this->initializerListQuad.rowsize() + j),
                      this->initializerListQuad(i, j));
}

/* TYPED_TEST(TestSparseMatrix, ValueConstIteratorTest) { */
/*     for (const auto& row : this->sizedConstructed) */
/*         for (const auto& element : row) EXPECT_EQ(TypeParam(0), element); */

/*     for (const auto& row : this->sizedConstructedZero) */
/*         for (const auto& element : row) EXPECT_EQ(TypeParam(0), element); */

/*     for (const auto& row : this->sizedConstructedValue) */
/*         for (const auto& element : row) */
/*             EXPECT_EQ(TypeParam(this->value), element); */
/* } */

/* TYPED_TEST(TestSparseMatrix, ValueIteratorTest) { */
/*     for (auto& row : this->sizedConstructed) */
/*         for (auto& element : row) EXPECT_EQ(TypeParam(0), element); */

/*     for (auto& row : this->sizedConstructedZero) */
/*         for (auto& element : row) EXPECT_EQ(TypeParam(0), element); */

/*     for (auto& row : this->sizedConstructedValue) */
/*         for (auto& element : row) EXPECT_EQ(TypeParam(this->value), element);
 */
/* } */

}  // namespace
