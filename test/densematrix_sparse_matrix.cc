#include <complex>
#include <initializer_list>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the copy constructor of DenseMatrix from a
// sparsematrix

using hdnum::SparseMatrix, hdnum::Vector, std::initializer_list;

template <typename T>
class TestDenseMatrixCopyConstructor : public ::testing::Test {
public:
    SparseMatrix<T> A;

    TestDenseMatrixCopyConstructor() {
        auto builder = typename SparseMatrix<T>::builder(4, 5);

        builder.addEntry(3, 2, 11);
        builder.addEntry(0, 0, 10);
        builder.addEntry(1, 2, 11);
        builder.addEntry(0, 3, 12);
        builder.addEntry(2, 1, 16);
        builder.addEntry(1, 4, 13);
        builder.addEntry(3, 4, 13);

        A = builder.build();
    }
};

using TestTypes = ::testing::Types<int, double, float, std::complex<int>,
                                   std::complex<double>, std::complex<float>>;

TYPED_TEST_SUITE(TestDenseMatrixCopyConstructor, TestTypes);

TYPED_TEST(TestDenseMatrixCopyConstructor, CopyConstructorFromSparseMatrix) {
    hdnum::DenseMatrix<TypeParam> B(this->A);

    ASSERT_EQ(this->A.rowsize(), B.rowsize());
    ASSERT_EQ(this->A.colsize(), B.colsize());

    for (std::size_t i = 0; i < B.rowsize(); i++) {
        for (std::size_t j = 0; j < B.colsize(); j++) {
            ASSERT_EQ(this->A(i, j), B(i, j));
        }
    }
}

}  // namespace
