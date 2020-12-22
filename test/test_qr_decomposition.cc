#include "gtest/gtest.h"
#include <ctime>

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
    // check that there's an exception by runnig qr decomposition
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

    // A = Q*R
    hdnum::DenseMatrix<double> QR(Q*R);
    ASSERT_EQ(QR.rowsize(), A.rowsize());
    ASSERT_EQ(QR.colsize(), A.colsize());
    for (int i=0; i < QR.rowsize(); i++) {
        for (int j=0; j < QR.colsize(); j++) {
            ASSERT_NEAR(QR(i, j), A(i, j), 0.00000001);
        }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
