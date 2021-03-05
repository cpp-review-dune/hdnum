#include <complex>
#include <initializer_list>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the operators for the SparseMatrix class.

using hdnum::SparseMatrix, hdnum::Vector, std::initializer_list;

template <typename T>
class TestSparseMatrixOperators : public ::testing::Test {
public:
    const SparseMatrix<T> A;
    const Vector<T> x;

    using size_type = typename SparseMatrix<T>::size_type;

    constexpr static inline const initializer_list<initializer_list<T>>
        initializerListQuad = {{1, 3, 4}, {-2, -4, 2}, {-1, 2, 2}};

    TestSparseMatrixOperators()
        : A(typename SparseMatrix<T>::builder(initializerListQuad).build()),
          x({-18, -20, 15}) {}
};

using TestTypes = ::testing::Types<int, double, float, std::complex<int>,
                                   std::complex<double>, std::complex<float>>;

TYPED_TEST_SUITE(TestSparseMatrixOperators, TestTypes);

TYPED_TEST(TestSparseMatrixOperators, ComparisonOperators) {
    SparseMatrix<TypeParam> B(this->A);
    EXPECT_TRUE(this->A == B);
    EXPECT_FALSE(this->A != B);

    auto A_t = this->A.transpose();

    EXPECT_FALSE(this->A == A_t);
    EXPECT_TRUE(this->A != A_t);
    EXPECT_TRUE(this->A == A_t.transpose());
}

TYPED_TEST(TestSparseMatrixOperators, AccessOperators) {
    SparseMatrix<TypeParam> B(this->A);
    EXPECT_EQ(B(0, 0), TypeParam {1});
    EXPECT_EQ(B(0, 1), TypeParam {3});
    EXPECT_EQ(B(0, 2), TypeParam {4});

    // B(0, 1) = TypeParam {5};

    // EXPECT_EQ(B(0, 1), TypeParam {5});
}

TYPED_TEST(TestSparseMatrixOperators, MatrixVectorMultiplication) {
    using size_type = typename TestSparseMatrixOperators<TypeParam>::size_type;

    const auto expectedResult = Vector<TypeParam>({-18, 146, 8});

    auto y = this->A * this->x;
    Vector<TypeParam> y_mv(this->A.rowsize());

    this->A.mv(y_mv, this->x);

    for (auto i = size_type(0); i < expectedResult.size(); i++) {
        EXPECT_EQ(y.at(i), expectedResult.at(i));
        EXPECT_EQ(y_mv.at(i), expectedResult.at(i));
    }
}

}  // namespace
