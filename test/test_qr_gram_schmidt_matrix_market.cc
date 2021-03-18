#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the qr decomposition via gram-schmidt process for
// sparse matrices

using hdnum::DenseMatrix;

class TestQRGramSchmidtSparse : public ::testing::Test {
public:
    const std::string filename_small_full_rank =
        "../matrix_market_files/example_small_full_rank.mtx";

    const std::string filename_square_full_rank =
        "../matrix_market_files/example_square_full_rank.mtx";

    const double threshold = 0.000000001;
};

// Test qr_gram_schmidt and qr_gram_schmidt_simple for a square matrix
TEST_F(TestQRGramSchmidtSparse, TestQRSquareMatrix) {
    DenseMatrix<double> A;
    hdnum::readMatrixFromFileMatrixMarket(filename_square_full_rank, A);

    DenseMatrix<double> Q(A);
    DenseMatrix<double> R(hdnum::qr_gram_schmidt(Q));
    DenseMatrix<double> QR(Q * R);

    DenseMatrix<double> Q_simple(A);
    DenseMatrix<double> R_simple(hdnum::qr_gram_schmidt_simple(Q_simple));
    DenseMatrix<double> QR_simple(Q_simple * R_simple);

    ASSERT_EQ(A.colsize(), QR.colsize());
    ASSERT_EQ(A.rowsize(), QR.rowsize());
    ASSERT_EQ(A.colsize(), QR_simple.colsize());
    ASSERT_EQ(A.rowsize(), QR_simple.rowsize());

    for (int i = 0; i < QR.rowsize(); i++) {
        for (int j = 0; j < QR.colsize(); j++) {
            ASSERT_NEAR(QR(i, j), A(i, j), threshold);
        }
    }

    for (int i = 0; i < QR_simple.rowsize(); i++) {
        for (int j = 0; j < QR_simple.colsize(); j++) {
            ASSERT_NEAR(QR_simple(i, j), A(i, j), threshold);
        }
    }
}

// Test qr_gram_schmidt and qr_gram_schmidt_simple for a small matrix
TEST_F(TestQRGramSchmidtSparse, TestQRSmallMatrix) {
    DenseMatrix<double> A;

    hdnum::readMatrixFromFileMatrixMarket(filename_small_full_rank, A);

    DenseMatrix<double> Q(A);
    DenseMatrix<double> R(hdnum::qr_gram_schmidt(Q));
    DenseMatrix<double> QR(Q * R);

    DenseMatrix<double> Q_simple(A);
    DenseMatrix<double> R_simple(hdnum::qr_gram_schmidt_simple(Q_simple));
    DenseMatrix<double> QR_simple(Q_simple * R_simple);

    ASSERT_EQ(A.colsize(), QR.colsize());
    ASSERT_EQ(A.rowsize(), QR.rowsize());
    ASSERT_EQ(A.colsize(), QR_simple.colsize());
    ASSERT_EQ(A.rowsize(), QR_simple.rowsize());

    for (int i = 0; i < QR.rowsize(); i++) {
        for (int j = 0; j < QR.colsize(); j++) {
            ASSERT_NEAR(QR(i, j), A(i, j), threshold);
        }
    }

    for (int i = 0; i < QR_simple.rowsize(); i++) {
        for (int j = 0; j < QR_simple.colsize(); j++) {
            ASSERT_NEAR(QR_simple(i, j), A(i, j), threshold);
        }
    }
}
}  // namespace