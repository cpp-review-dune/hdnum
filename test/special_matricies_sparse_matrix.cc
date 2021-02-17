#include <complex>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test all the special matricies related to the
// SparseMatrix class.

using hdnum::SparseMatrix;

template <typename T>
class TestSparseMatrixSpecial : public ::testing::Test {
public:
    // This will be created by size
    const SparseMatrix<T> sizedIdentity;

    const SparseMatrix<T> templateMatrix;
    // this will be created so it matches the size of the template matrix
    const SparseMatrix<T> matchingIdentity;

    // this will overwrite a matrix to be an identity
    SparseMatrix<T> overwriteRefIdentity;

    using size_type = typename SparseMatrix<T>::size_type;

    constexpr static inline const size_type dimN = 5;
    constexpr static inline const size_type dimM = 5;

    TestSparseMatrixSpecial()
        : sizedIdentity(SparseMatrix<T>::identity(this->dimN)),
          templateMatrix(dimM, dimN),
          matchingIdentity(this->templateMatrix.matchingIdentity()),
          overwriteRefIdentity(dimM, dimN) {
        identity(overwriteRefIdentity);
    }
};

using TestTypes = ::testing::Types<int, double, float, std::complex<int>,
                                   std::complex<double>, std::complex<float>>;
TYPED_TEST_SUITE(TestSparseMatrixSpecial, TestTypes);

TYPED_TEST(TestSparseMatrixSpecial, SizeTest) {
    using T = TestSparseMatrixSpecial<TypeParam>;

    EXPECT_EQ(T::dimN, this->sizedIdentity.rowsize());
    EXPECT_EQ(T::dimN, this->sizedIdentity.colsize());

    EXPECT_EQ(T::dimN, this->matchingIdentity.rowsize());
    EXPECT_EQ(T::dimN, this->matchingIdentity.colsize());
    EXPECT_EQ(this->templateMatrix.colsize(), this->matchingIdentity.rowsize());
    EXPECT_EQ(this->templateMatrix.rowsize(), this->matchingIdentity.colsize());

    EXPECT_EQ(T::dimN, this->overwriteRefIdentity.rowsize());
    EXPECT_EQ(T::dimN, this->overwriteRefIdentity.colsize());
}

TYPED_TEST(TestSparseMatrixSpecial, ValueIndexTest) {
    using size_type = typename SparseMatrix<TypeParam>::size_type;

    for (auto i = size_type(0); this->sizedIdentity.rowsize(); i++)
        for (auto j = size_type(0); this->sizedIdentity.colsize(); j++)
            if (i != j)
                EXPECT_EQ(TypeParam(0), this->sizedIdentity(i, j));
            else
                EXPECT_EQ(TypeParam(1), this->sizedIdentity(i, j));

    for (auto i = size_type(0); this->matchingIdentity.rowsize(); i++)
        for (auto j = size_type(0); this->matchingIdentity.colsize(); j++)
            if (i != j)
                EXPECT_EQ(TypeParam(0), this->matchingIdentity(i, j));
            else
                EXPECT_EQ(TypeParam(1), this->matchingIdentity(i, j));

    for (auto i = size_type(0); this->overwriteRefIdentity.rowsize(); i++)
        for (auto j = size_type(0); this->overwriteRefIdentity.colsize(); j++)
            if (i != j)
                EXPECT_EQ(TypeParam(0), this->overwriteRefIdentity(i, j));
            else
                EXPECT_EQ(TypeParam(1), this->overwriteRefIdentity(i, j));
}

TYPED_TEST(TestSparseMatrixSpecial, VerifyMult) {
    using size_type = typename SparseMatrix<TypeParam>::size_type;

    SparseMatrix<TypeParam> testMatrix(this->dimM, this->dimN);
    for (auto i = size_type(0); testMatrix.rowsize(); i++)
        for (auto j = size_type(0); testMatrix.colsize(); j++)
            testMatrix.get(i, j) = TypeParam(i + j);
}

}  // namespace
