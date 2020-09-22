#include <complex>
#include <typeinfo>
#include <initializer_list>
#include <type_traits>

#include "hdnum.hh"
#include "gtest/gtest.h"

namespace {
// here we test the past-conditions of the qr decomposition

    template<typename T>
    class TestQRDecomposition : public ::testing::Test {
    public:
        const hdnum::DenseMatrix<T> A;
        hdnum::DenseMatrix<T> Q;
        hdnum::DenseMatrix<T> R;
        hdnum::DenseMatrix<T> QR;

        constexpr static inline const std::initializer_list< std::initializer_list<T> >
            random_initializer_list = {
                {-1, 1, 4, -1},
                {3, 8, 1, -4},
                {7, 3, -1, 2},
                {2, -4, -1, 6}
        };

        TestQRDecomposition() : A(random_initializer_list), Q(4, 4), R(4, 4), QR(4, 4) {
            hdnum::qr_decomposition_gram_schmidt(A, Q, R);
            QR = Q*R;
        }
    };

    using TestTypes = ::testing::Types< float, double, std::complex<double>>;

    TYPED_TEST_SUITE(TestQRDecomposition, TestTypes);

    TYPED_TEST(TestQRDecomposition, TestValues) {
        for (int i=0; i < this->A.rowsize(); i++) {
            for (int j=0; j < this->A.colsize(); j++) {
                if constexpr(std::is_same_v<double, TypeParam>)
                    EXPECT_DOUBLE_EQ(this->A[i][j], this->QR[i][j]);
                else if constexpr(std::is_same_v<float, TypeParam>)
                    EXPECT_FLOAT_EQ(this->A[i][j], this->QR[i][j]);
                else if constexpr(std::is_same_v<std::complex<float>, TypeParam>) {
                    EXPECT_FLOAT_EQ(this->A[i][j].real(), this->QR[i][j].real());
                    EXPECT_FLOAT_EQ(this->A[i][j].imag(), this->QR[i][j].imag());
                }
                else if constexpr(std::is_same_v<std::complex<double>, TypeParam>) {
                    EXPECT_DOUBLE_EQ(this->A[i][j].real(), this->QR[i][j].real());
                    EXPECT_DOUBLE_EQ(this->A[i][j].imag(), this->QR[i][j].imag());
                }
                else
                    EXPECT_EQ(this->A[i][j], this->QR[i][j]);
            }
        }
    }

} // namespace
