#include "gtest/gtest.h"
#include <ctime>
#include <gmpxx.h>

#include "../hdnum.hh"


/*   What has to be tested? 
 *
 *   Preconditions:
 *    - A is a small matrix (A ∈ K^(mxn), m > n), otherwise exception
 *
 *   Postconditions:
 *    - Q ∈ K^(mxn)
 *    - Q is a orthogonal matrix: 
 *        -> ||q1|| = 1, for all q1, ... , qn
 *        -> Q^(T) * Q = I
 *    - R ∈ K^(nxn)
 *    - R is an upper triangular matrix
 *    - Q*R = A ∈ K^(mxn)
 */


// check preconditions
TEST(TestQRDecomposition, TestPreconditions)
{
    // create a square matrix
    hdnum::DenseMatrix<double> A(10, 10);
    // check that there's an exception by running qr decomposition
    EXPECT_THROW({
        try
        {
            hdnum::qr_decomposition_gram_schmidt(A);
        }
        catch(const hdnum::ErrorException& e)
        {
            // check that the exception message is right
            std::size_t found = e.what().find("The Matrix is not a small matrix!");
            ASSERT_TRUE(found != e.what().npos);
            throw;
        }
    }, hdnum::ErrorException);

    // repeat with a matrix n > m
    hdnum::DenseMatrix<double> B(6, 10);
    // check that there's an exception by runnig qr decomposition
    EXPECT_THROW({
        try
        {
            hdnum::qr_decomposition_gram_schmidt(B);
        }
        catch(const hdnum::ErrorException& e)
        {
            // check that the exception message is right
            std::size_t found = e.what().find("The Matrix is not a small matrix!");
            ASSERT_TRUE(found != e.what().npos);
            throw;
        }
    }, hdnum::ErrorException);
    
}


// check postconditions
TEST(TestQRDecomposition, TestPostconditions)
{
    // create a random size small matrix
    int m;
    int n;
    do {
        srand(time(NULL));
        m = (rand() % 18) + 2;    // [2, ... ,20]
        n = (rand() % 18) + 2;    // [2, ... ,20]
    } while (m <= n);
    hdnum::DenseMatrix<double> Q(m, n);

    // fill it with random elements
    for (int i=0; i < Q.rowsize(); i++) {
        for (int j=0; j < Q.colsize(); j++) {
            int x = (rand() % 200) - 100;    // [-100, ... ,100]
            Q(i, j) = x;
        }
    }

    // save it before overwritting
    hdnum::DenseMatrix<double> A(Q);

    // run qr decomposition and save R
    hdnum::DenseMatrix<double> R(hdnum::qr_decomposition_gram_schmidt(Q));

    // declare a threshold for the asserts
    double threshold = 0.00000001;

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q.colsize(), A.colsize());
    ASSERT_EQ(Q.rowsize(), A.rowsize());

    // check that norm(qi) = 1
    for (int i=0; i < Q.colsize(); i++) {
        double norm(0.0);

        for (int j=0; j < Q.rowsize(); j++) {
            norm += Q(j, i) * Q(j, i);
        }
        norm = sqrt(norm);
        ASSERT_NEAR(norm, 1.0, threshold);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I(Q.transpose()*Q);
    for (int i=0; i < I.rowsize(); i++) {
        for (int j=0; j < I.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I(i, j), 0.0, threshold);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R.rowsize(), A.colsize());
    ASSERT_EQ(R.colsize(), A.colsize());

    // R is an upper triangular matrix
    for (int i=0; i < R.colsize(); i++) {
        for (int j=i+1; j < R.rowsize(); j++) {
            ASSERT_NEAR(R(j, i), 0.0, threshold);
        }
    }

    // A = Q*R
    hdnum::DenseMatrix<double> QR(Q*R);
    ASSERT_EQ(QR.rowsize(), A.rowsize());
    ASSERT_EQ(QR.colsize(), A.colsize());
    for (int i=0; i < QR.rowsize(); i++) {
        for (int j=0; j < QR.colsize(); j++) {
            ASSERT_NEAR(QR(i, j), A(i, j), threshold);
        }
    }

    // now repeat with GNU multiple precision types
    mpf_set_default_prec(2048);

    // declare a new more precise threshold
    mpf_class threshold_gmp = 0.00000000000000000001;

    // create A, Q and R with gmp datatypes
    hdnum::DenseMatrix<mpf_class> A_gmp(m, n);
    for (int r=0; r < A_gmp.rowsize(); r++) {
        for (int c=0; c < A_gmp.colsize(); c++) {
            A_gmp(r, c) = A(r, c);
        }
    }

    hdnum::DenseMatrix<mpf_class> Q_gmp(A_gmp);
    hdnum::DenseMatrix<mpf_class> R_gmp(hdnum::qr_decomposition_gram_schmidt(Q_gmp));

    // check that norm(qi) = 1
    for (int i=0; i < Q_gmp.colsize(); i++) {
        mpf_class norm = 0;
        mpf_class expected_norm = 1;

        for (int j=0; j < Q_gmp.rowsize(); j++) {
            norm += Q_gmp(j, i) * Q_gmp(j, i);
        }
        norm = sqrt(norm);
        ASSERT_TRUE( (norm > expected_norm-threshold_gmp) && (norm < expected_norm+threshold_gmp));     // 1-threshold < norm < 1+threshold
    }

    // check that Q_gmp^T * Q_gmp = I
    hdnum::DenseMatrix<mpf_class> I_gmp(Q_gmp.transpose()*Q_gmp);
    for (int i=0; i < I_gmp.rowsize(); i++) {
        for (int j=0; j < I_gmp.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_TRUE((I_gmp(i, i) > 1-threshold_gmp) && (I_gmp(i, i) < 1+threshold_gmp));         // 1-threshold < I(i, i) < 1+threshold
                continue;
            }
            // other elements
            ASSERT_TRUE((I_gmp(i, j) > 0-threshold_gmp) && (I_gmp(i, j) < 0+threshold_gmp));             // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // R_gmp ∈ K^(nxn)
    ASSERT_EQ(R_gmp.rowsize(), A_gmp.colsize());
    ASSERT_EQ(R_gmp.colsize(), A_gmp.colsize());

    // R_gmp is an upper triangular matrix
    for (int i=0; i < R_gmp.colsize(); i++) {
        for (int j=i+1; j < R_gmp.rowsize(); j++) {
            ASSERT_TRUE((R_gmp(j, i) > 0-threshold_gmp) && (R_gmp(j, i) < 0+threshold_gmp));             // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // A_gmp = Q_gmp*R_gmp
    hdnum::DenseMatrix<mpf_class> QR_gmp(Q_gmp*R_gmp);
    ASSERT_EQ(QR_gmp.rowsize(), A_gmp.rowsize());
    ASSERT_EQ(QR_gmp.colsize(), A_gmp.colsize());
    for (int i=0; i < QR_gmp.rowsize(); i++) {
        for (int j=0; j < QR_gmp.colsize(); j++) {
            ASSERT_TRUE((QR_gmp(i, j) > A_gmp(i, j)-threshold_gmp) && (QR_gmp(i, j) < A_gmp(i, j)+threshold_gmp));     // A_gmp(i, j)-threshold < QR_gmp(i, j) < A_gmp(i, j)+threshold
        }
    }
}

// check preconditions
TEST(TestQRDecompositionPivoting, TestPreconditions)
{
    // create a square matrix
    hdnum::DenseMatrix<double> A(10, 10);

    // check that there's an exception by runnig qr decomposition
    EXPECT_THROW({
        try
        {
            hdnum::qr_decomposition_gram_schmidt_pivoting(A);
        }
        catch(const hdnum::ErrorException& e)
        {
            // check that the exception message is right
            std::size_t found = e.what().find("The Matrix is not a small matrix!");
            ASSERT_TRUE(found != e.what().npos);
            throw;
        }
    }, hdnum::ErrorException);

    // repeat with a matrix n > m
    hdnum::DenseMatrix<double> B(6, 10);
    // check that there's an exception by runnig qr decomposition
    EXPECT_THROW({
        try
        {
            hdnum::qr_decomposition_gram_schmidt_pivoting(B);
        }
        catch(const hdnum::ErrorException& e)
        {
            // check that the exception message is right
            std::size_t found = e.what().find("The Matrix is not a small matrix!");
            ASSERT_TRUE(found != e.what().npos);
            throw;
        }
    }, hdnum::ErrorException);
    
}

// check postconditions
TEST(TestQRDecompositionPivoting, TestPostconditions)
{
    // create a random size small matrix
    int m;
    int n;
    do {
        srand(time(NULL));
        m = (rand() % 18) + 2;    // [2, ... ,20]
        n = (rand() % 18) + 2;    // [2, ... ,20]
    } while (m <= n);
    hdnum::DenseMatrix<double> Q(m, n);

    // fill it with random elements
    for (int i=0; i < Q.rowsize(); i++) {
        for (int j=0; j < Q.colsize(); j++) {
            int x = (rand() % 200) - 100;    // [-100, ... ,100]
            Q(i, j) = x;
        }
    }

    // save it before overwritting
    hdnum::DenseMatrix<double> A(Q);

    // run qr decomposition and save R
    hdnum::DenseMatrix<double> R(hdnum::qr_decomposition_gram_schmidt_pivoting(Q));

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q.colsize(), A.colsize());
    ASSERT_EQ(Q.rowsize(), A.rowsize());

    // check that norm(qi) = 1
    for (int i=0; i < Q.colsize(); i++) {
        double norm(0.0);

        for (int j=0; j < Q.rowsize(); j++) {
            norm += Q(j, i) * Q(j, i);
        }
        norm = sqrt(norm);
        ASSERT_NEAR(norm, 1.0, 0.00000001);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I(Q.transpose()*Q);
    for (int i=0; i < I.rowsize(); i++) {
        for (int j=0; j < I.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I(i, j), 1.0, 0.00000001);
                continue;
            }
            // other elements
            ASSERT_NEAR(I(i, j), 0.0, 0.00000001);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R.rowsize(), A.colsize());
    ASSERT_EQ(R.colsize(), A.colsize());

    // R is an upper triangular matrix
    for (int i=0; i < R.colsize(); i++) {
        for (int j=i+1; j < R.rowsize(); j++) {
            ASSERT_NEAR(R(j, i), 0.0, 0.00000001);
        }
    }

    // now repeat with GNU multiple precision types
    mpf_set_default_prec(2048);

    // declare a new more precise threshold
    mpf_class threshold_gmp = 0.00000000000000000001;

    // create A, Q and R with gmp datatypes
    hdnum::DenseMatrix<mpf_class> A_gmp(m, n);
    for (int r=0; r < A_gmp.rowsize(); r++) {
        for (int c=0; c < A_gmp.colsize(); c++) {
            A_gmp(r, c) = A(r, c);
        }
    }

    hdnum::DenseMatrix<mpf_class> Q_gmp(A_gmp);
    hdnum::DenseMatrix<mpf_class> R_gmp(hdnum::qr_decomposition_gram_schmidt_pivoting(Q_gmp));

    // check that norm(qi) = 1
    for (int i=0; i < Q_gmp.colsize(); i++) {
        mpf_class norm = 0;
        mpf_class expected_norm = 1;

        for (int j=0; j < Q_gmp.rowsize(); j++) {
            norm += Q_gmp(j, i) * Q_gmp(j, i);
        }
        norm = sqrt(norm);
        ASSERT_TRUE( (norm > expected_norm-threshold_gmp) && (norm < expected_norm+threshold_gmp));     // 1-threshold < norm < 1+threshold
    }

    // check that Q_gmp^T * Q_gmp = I
    hdnum::DenseMatrix<mpf_class> I_gmp(Q_gmp.transpose()*Q_gmp);
    for (int i=0; i < I_gmp.rowsize(); i++) {
        for (int j=0; j < I_gmp.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_TRUE((I_gmp(i, i) > 1-threshold_gmp) && (I_gmp(i, i) < 1+threshold_gmp));         // 1-threshold < I(i, i) < 1+threshold
                continue;
            }
            // other elements
            ASSERT_TRUE((I_gmp(i, j) > 0-threshold_gmp) && (I_gmp(i, j) < 0+threshold_gmp));             // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // R_gmp ∈ K^(nxn)
    ASSERT_EQ(R_gmp.rowsize(), A_gmp.colsize());
    ASSERT_EQ(R_gmp.colsize(), A_gmp.colsize());

    // R_gmp is an upper triangular matrix
    for (int i=0; i < R_gmp.colsize(); i++) {
        for (int j=i+1; j < R_gmp.rowsize(); j++) {
            ASSERT_TRUE((R_gmp(j, i) > 0-threshold_gmp) && (R_gmp(j, i) < 0+threshold_gmp));             // 0-threshold < I(i, j) < 0+threshold
        }
    }
}

// check rank reveal
TEST(TestQRDecompositionPivoting, TestRankReveal)
{
    // create a Matrix A with rank(A) < number of columns (a0 and a2 are linear dependant)
    hdnum::DenseMatrix<double> A(5, 4);

    A(0, 0) = 1;     A(0, 1) = 3;    A(0, 2) = 5;     A(0, 3) = 6;
    A(1, 0) = 3;     A(1, 1) = 3;    A(1, 2) = 15;    A(1, 3) = 5;
    A(2, 0) = 2;     A(2, 1) = 7;    A(2, 2) = 10;    A(2, 3) = 6;
    A(3, 0) = 5;     A(3, 1) = 3;    A(3, 2) = 25;    A(3, 3) = 0;
    A(4, 0) = 11;    A(4, 1) = 6;    A(4, 2) = 55;    A(4, 3) = 6;

    // check that there's an exception by runnig qr decomposition
    EXPECT_THROW({
        try
        {
            hdnum::qr_decomposition_gram_schmidt_pivoting(A);
        }
        catch(const hdnum::ErrorException& e)
        {
            // check that the exception message is right
            std::size_t found = e.what().find("The Matrix has not full rank!");
            ASSERT_TRUE(found != e.what().npos);
            throw;
        }
    }, hdnum::ErrorException);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
