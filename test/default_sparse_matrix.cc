#include <complex>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the empty SparseMatrix class.

using hdnum::SparseMatrix;

template <typename T>
class TestSparseMatrix : public ::testing::Test {
public:
    const hdnum::SparseMatrix<T> s1;
    const hdnum::SparseMatrix<T> s2;
    TestSparseMatrix() : s1(), s2(s1) {}
};

using TestTypes = ::testing::Types<int, double, float, std::complex<int>,
                                   std::complex<double>, std::complex<float>>;
TYPED_TEST_SUITE(TestSparseMatrix, TestTypes);

TYPED_TEST(TestSparseMatrix, SizeTest) {
    using size_type = typename SparseMatrix<TypeParam>::size_type;

    // default-constructed matrix should have 0x0
    EXPECT_EQ(size_type(0), this->s1.rowsize());
    EXPECT_EQ(size_type(0), this->s1.colsize());

    // copy-constructed matrix should also have same dimensions
    EXPECT_EQ(size_type(0), this->s2.rowsize());
    EXPECT_EQ(size_type(0), this->s2.colsize());
}

}  // namespace
