#include <gmpxx.h>

#include <ctime>
#include <iostream>

#include "../hdnum.hh"
#include "gtest/gtest.h"

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

using size_type = hdnum::DenseMatrix<int>::size_type;

// check postconditions
TEST(TestQRDecomposition, TestPostconditionsSmallMatrix) {
    // create a random size small matrix
    int m;
    int n;
    do {
        srand(time(NULL));
        m = (rand() % 18) + 2;  // [2, ... ,20]
        n = (rand() % 18) + 2;  // [2, ... ,20]
    } while (m == n);
    if (n > m) {
        std::swap(m, n);
    }
    hdnum::DenseMatrix<double> Q1(m, n);

    // fill it with random elements
    for (int i = 0; i < Q1.rowsize(); i++) {
        for (int j = 0; j < Q1.colsize(); j++) {
            int x = (rand() % 200) - 100;  // [-100, ... ,100]
            Q1(i, j) = x;
        }
    }

    // save it before overwritting
    hdnum::DenseMatrix<double> A(Q1);

    // copy it for another qr function
    hdnum::DenseMatrix<double> Q2(Q1);

    // run qr decomposition and save R
    hdnum::DenseMatrix<double> R1(hdnum::qr_gram_schmidt_simple(Q1));
    hdnum::DenseMatrix<double> R2(hdnum::qr_gram_schmidt(Q2));

    // declare a threshold for the asserts
    double threshold = 0.00000001;

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q1.colsize(), A.colsize());
    ASSERT_EQ(Q1.rowsize(), A.rowsize());
    ASSERT_EQ(Q2.colsize(), A.colsize());
    ASSERT_EQ(Q2.rowsize(), A.rowsize());

    // check that norm(qi) = 1
    for (int i = 0; i < Q1.colsize(); i++) {
        double normQ1(0.0);
        double normQ2(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ1 += Q1(j, i) * Q1(j, i);
            normQ2 += Q2(j, i) * Q2(j, i);
        }
        normQ1 = sqrt(normQ1);
        normQ2 = sqrt(normQ2);
        ASSERT_NEAR(normQ1, 1.0, threshold);
        ASSERT_NEAR(normQ2, 1.0, threshold);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I1(Q1.transpose() * Q1);
    hdnum::DenseMatrix<double> I2(Q2.transpose() * Q2);
    for (int i = 0; i < I1.rowsize(); i++) {
        for (int j = 0; j < I1.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I1(i, j), 1.0, threshold);
                ASSERT_NEAR(I2(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I1(i, j), 0.0, threshold);
            ASSERT_NEAR(I2(i, j), 0.0, threshold);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R1.rowsize(), A.colsize());
    ASSERT_EQ(R1.colsize(), A.colsize());
    ASSERT_EQ(R2.rowsize(), A.colsize());
    ASSERT_EQ(R2.colsize(), A.colsize());

    // R is an upper triangular matrix
    for (int i = 0; i < R1.colsize(); i++) {
        for (int j = i + 1; j < R1.rowsize(); j++) {
            ASSERT_NEAR(R1(j, i), 0.0, threshold);
            ASSERT_NEAR(R2(j, i), 0.0, threshold);
        }
    }

    // A = Q*R
    hdnum::DenseMatrix<double> QR1(Q1 * R1);
    hdnum::DenseMatrix<double> QR2(Q2 * R2);
    ASSERT_EQ(QR1.rowsize(), A.rowsize());
    ASSERT_EQ(QR1.colsize(), A.colsize());
    ASSERT_EQ(QR2.rowsize(), A.rowsize());
    ASSERT_EQ(QR2.colsize(), A.colsize());
    for (int i = 0; i < QR1.rowsize(); i++) {
        for (int j = 0; j < QR1.colsize(); j++) {
            ASSERT_NEAR(QR1(i, j), A(i, j), threshold);
            ASSERT_NEAR(QR2(i, j), A(i, j), threshold);
        }
    }

    // now repeat with GNU multiple precision types
    mpf_set_default_prec(2048);

    // declare a new more precise threshold
    mpf_class threshold_gmp = 0.00000000000000000001;

    // create A, Q and R with gmp datatypes
    hdnum::DenseMatrix<mpf_class> A_gmp(m, n);
    for (int r = 0; r < A_gmp.rowsize(); r++) {
        for (int c = 0; c < A_gmp.colsize(); c++) {
            A_gmp(r, c) = A(r, c);
        }
    }

    hdnum::DenseMatrix<mpf_class> Q1_gmp(A_gmp);
    hdnum::DenseMatrix<mpf_class> Q2_gmp(A_gmp);

    hdnum::DenseMatrix<mpf_class> R1_gmp(
        hdnum::qr_gram_schmidt(Q1_gmp));
    hdnum::DenseMatrix<mpf_class> R2_gmp(
        hdnum::qr_gram_schmidt(Q2_gmp));

    // check that norm(qi) = 1
    for (int i = 0; i < Q1_gmp.colsize(); i++) {
        mpf_class normQ1 = 0;
        mpf_class normQ2 = 0;
        mpf_class expected_norm = 1;

        for (int j = 0; j < Q1_gmp.rowsize(); j++) {
            normQ1 += Q1_gmp(j, i) * Q1_gmp(j, i);
            normQ2 += Q2_gmp(j, i) * Q2_gmp(j, i);
        }
        normQ1 = sqrt(normQ1);
        normQ2 = sqrt(normQ2);
        ASSERT_TRUE(
            (normQ1 > expected_norm - threshold_gmp) &&
            (normQ1 < expected_norm +
                        threshold_gmp));  // 1-threshold < norm < 1+threshold
        ASSERT_TRUE(
            (normQ2 > expected_norm - threshold_gmp) &&
            (normQ2 < expected_norm +
                        threshold_gmp));  // 1-threshold < norm < 1+threshold
    }

    // check that Q_gmp^T * Q_gmp = I
    hdnum::DenseMatrix<mpf_class> I1_gmp(Q1_gmp.transpose() * Q1_gmp);
    hdnum::DenseMatrix<mpf_class> I2_gmp(Q2_gmp.transpose() * Q2_gmp);
    for (int i = 0; i < I1_gmp.rowsize(); i++) {
        for (int j = 0; j < I1_gmp.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_TRUE(
                    (I1_gmp(i, i) > 1 - threshold_gmp) &&
                    (I1_gmp(i, i) < 1 + threshold_gmp));
                ASSERT_TRUE(
                    (I2_gmp(i, i) > 1 - threshold_gmp) &&
                    (I2_gmp(i, i) < 1 + threshold_gmp));
                continue;
            }
            // other elements
            ASSERT_TRUE(
                (I1_gmp(i, j) > 0 - threshold_gmp) &&
                (I1_gmp(i, j) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
            ASSERT_TRUE(
                (I2_gmp(i, j) > 0 - threshold_gmp) &&
                (I2_gmp(i, j) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // R_gmp ∈ K^(nxn)
    ASSERT_EQ(R1_gmp.rowsize(), A_gmp.colsize());
    ASSERT_EQ(R1_gmp.colsize(), A_gmp.colsize());
    ASSERT_EQ(R2_gmp.rowsize(), A_gmp.colsize());
    ASSERT_EQ(R2_gmp.colsize(), A_gmp.colsize());

    // R_gmp is an upper triangular matrix
    for (int i = 0; i < R1_gmp.colsize(); i++) {
        for (int j = i + 1; j < R1_gmp.rowsize(); j++) {
            ASSERT_TRUE(
                (R1_gmp(j, i) > 0 - threshold_gmp) &&
                (R1_gmp(j, i) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
            ASSERT_TRUE(
                (R2_gmp(j, i) > 0 - threshold_gmp) &&
                (R2_gmp(j, i) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // A_gmp = Q_gmp*R_gmp
    hdnum::DenseMatrix<mpf_class> QR1_gmp(Q1_gmp * R1_gmp);
    hdnum::DenseMatrix<mpf_class> QR2_gmp(Q2_gmp * R2_gmp);
    ASSERT_EQ(QR1_gmp.rowsize(), A_gmp.rowsize());
    ASSERT_EQ(QR1_gmp.colsize(), A_gmp.colsize());
    ASSERT_EQ(QR2_gmp.rowsize(), A_gmp.rowsize());
    ASSERT_EQ(QR2_gmp.colsize(), A_gmp.colsize());
    for (int i = 0; i < QR1_gmp.rowsize(); i++) {
        for (int j = 0; j < QR1_gmp.colsize(); j++) {
            ASSERT_TRUE(
                (QR1_gmp(i, j) > A_gmp(i, j) - threshold_gmp) &&
                (QR1_gmp(i, j) <
                 A_gmp(i, j) +
                     threshold_gmp));
            ASSERT_TRUE(
                (QR2_gmp(i, j) > A_gmp(i, j) - threshold_gmp) &&
                (QR2_gmp(i, j) <
                 A_gmp(i, j) +
                     threshold_gmp));
        }
    }
}

// check postconditions
TEST(TestQRDecomposition, TestPostconditionsSquareMatrix) {
    // create a random size small matrix
    int m;
    srand(time(NULL));
    m = (rand() % 18) + 2;  // [2, ... ,20]
    hdnum::DenseMatrix<double> Q1(m, m);

    // fill it with random elements
    for (int i = 0; i < Q1.rowsize(); i++) {
        for (int j = 0; j < Q1.colsize(); j++) {
            int x = (rand() % 200) - 100;  // [-100, ... ,100]
            Q1(i, j) = x;
        }
    }

    // save it before overwritting
    hdnum::DenseMatrix<double> A(Q1);

    // create another Q for the second qr method
    hdnum::DenseMatrix<double> Q2(Q1);

    // run qr decomposition and save R
    hdnum::DenseMatrix<double> R1(hdnum::qr_gram_schmidt_simple(Q1));
    hdnum::DenseMatrix<double> R2(hdnum::qr_gram_schmidt(Q2));

    // declare a threshold for the asserts
    double threshold = 0.00000001;

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q1.colsize(), A.colsize());
    ASSERT_EQ(Q1.rowsize(), A.rowsize());
    ASSERT_EQ(Q2.colsize(), A.colsize());
    ASSERT_EQ(Q2.rowsize(), A.rowsize());

    // check that norm(qi) = 1
    for (int i = 0; i < Q1.colsize(); i++) {
        double normQ1(0.0);
        double normQ2(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ1 += Q1(j, i) * Q1(j, i);
            normQ2 += Q2(j, i) * Q2(j, i);
        }
        normQ1 = sqrt(normQ1);
        normQ2 = sqrt(normQ2);
        ASSERT_NEAR(normQ1, 1.0, threshold);
        ASSERT_NEAR(normQ2, 1.0, threshold);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I1(Q1.transpose() * Q1);
    hdnum::DenseMatrix<double> I2(Q2.transpose() * Q2);
    for (int i = 0; i < I1.rowsize(); i++) {
        for (int j = 0; j < I1.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I1(i, j), 1.0, threshold);
                ASSERT_NEAR(I2(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I1(i, j), 0.0, threshold);
            ASSERT_NEAR(I2(i, j), 0.0, threshold);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R1.rowsize(), A.colsize());
    ASSERT_EQ(R1.colsize(), A.colsize());
    ASSERT_EQ(R2.rowsize(), A.colsize());
    ASSERT_EQ(R2.colsize(), A.colsize());

    // R is an upper triangular matrix
    for (int i = 0; i < R1.colsize(); i++) {
        for (int j = i + 1; j < R1.rowsize(); j++) {
            ASSERT_NEAR(R1(j, i), 0.0, threshold);
            ASSERT_NEAR(R2(j, i), 0.0, threshold);
        }
    }

    // A = Q*R
    hdnum::DenseMatrix<double> QR1(Q1 * R1);
    hdnum::DenseMatrix<double> QR2(Q2 * R2);
    ASSERT_EQ(QR1.rowsize(), A.rowsize());
    ASSERT_EQ(QR1.colsize(), A.colsize());
    ASSERT_EQ(QR2.rowsize(), A.rowsize());
    ASSERT_EQ(QR2.colsize(), A.colsize());
    for (int i = 0; i < QR1.rowsize(); i++) {
        for (int j = 0; j < QR1.colsize(); j++) {
            ASSERT_NEAR(QR1(i, j), A(i, j), threshold);
            ASSERT_NEAR(QR2(i, j), A(i, j), threshold);
        }
    }

    // now repeat with GNU multiple precision types
    mpf_set_default_prec(2048);

    // declare a new more precise threshold
    mpf_class threshold_gmp = 0.00000000000000000001;

    // create A, Q and R with gmp datatypes
    hdnum::DenseMatrix<mpf_class> A_gmp(m, m);
    for (int r = 0; r < A_gmp.rowsize(); r++) {
        for (int c = 0; c < A_gmp.colsize(); c++) {
            A_gmp(r, c) = A(r, c);
        }
    }

    hdnum::DenseMatrix<mpf_class> Q1_gmp(A_gmp);
    hdnum::DenseMatrix<mpf_class> Q2_gmp(A_gmp);

    hdnum::DenseMatrix<mpf_class> R1_gmp(
        hdnum::qr_gram_schmidt(Q1_gmp));
    hdnum::DenseMatrix<mpf_class> R2_gmp(
        hdnum::qr_gram_schmidt(Q2_gmp));

    // check that norm(qi) = 1
    for (int i = 0; i < Q1_gmp.colsize(); i++) {
        mpf_class normQ1 = 0;
        mpf_class normQ2 = 0;
        mpf_class expected_norm = 1;

        for (int j = 0; j < Q1_gmp.rowsize(); j++) {
            normQ1 += Q1_gmp(j, i) * Q1_gmp(j, i);
            normQ2 += Q2_gmp(j, i) * Q2_gmp(j, i);
        }
        normQ1 = sqrt(normQ1);
        normQ2 = sqrt(normQ2);
        ASSERT_TRUE(
            (normQ1 > expected_norm - threshold_gmp) &&
            (normQ1 < expected_norm +
                        threshold_gmp));  // 1-threshold < norm < 1+threshold
        ASSERT_TRUE(
            (normQ2 > expected_norm - threshold_gmp) &&
            (normQ2 < expected_norm +
                        threshold_gmp));  // 1-threshold < norm < 1+threshold
    }

    // check that Q_gmp^T * Q_gmp = I
    hdnum::DenseMatrix<mpf_class> I1_gmp(Q1_gmp.transpose() * Q1_gmp);
    hdnum::DenseMatrix<mpf_class> I2_gmp(Q2_gmp.transpose() * Q2_gmp);
    for (int i = 0; i < I1_gmp.rowsize(); i++) {
        for (int j = 0; j < I1_gmp.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_TRUE(
                    (I1_gmp(i, i) > 1 - threshold_gmp) &&
                    (I1_gmp(i, i) < 1 + threshold_gmp));
                ASSERT_TRUE(
                    (I2_gmp(i, i) > 1 - threshold_gmp) &&
                    (I2_gmp(i, i) < 1 + threshold_gmp));
                continue;
            }
            // other elements
            ASSERT_TRUE(
                (I1_gmp(i, j) > 0 - threshold_gmp) &&
                (I1_gmp(i, j) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
            ASSERT_TRUE(
                (I2_gmp(i, j) > 0 - threshold_gmp) &&
                (I2_gmp(i, j) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // R_gmp ∈ K^(nxn)
    ASSERT_EQ(R1_gmp.rowsize(), A_gmp.colsize());
    ASSERT_EQ(R1_gmp.colsize(), A_gmp.colsize());
    ASSERT_EQ(R2_gmp.rowsize(), A_gmp.colsize());
    ASSERT_EQ(R2_gmp.colsize(), A_gmp.colsize());

    // R_gmp is an upper triangular matrix
    for (int i = 0; i < R1_gmp.colsize(); i++) {
        for (int j = i + 1; j < R1_gmp.rowsize(); j++) {
            ASSERT_TRUE(
                (R1_gmp(j, i) > 0 - threshold_gmp) &&
                (R1_gmp(j, i) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
            ASSERT_TRUE(
                (R2_gmp(j, i) > 0 - threshold_gmp) &&
                (R2_gmp(j, i) <
                 0 + threshold_gmp));  // 0-threshold < I(i, j) < 0+threshold
        }
    }

    // A_gmp = Q_gmp*R_gmp
    hdnum::DenseMatrix<mpf_class> QR1_gmp(Q1_gmp * R1_gmp);
    hdnum::DenseMatrix<mpf_class> QR2_gmp(Q2_gmp * R2_gmp);
    ASSERT_EQ(QR1_gmp.rowsize(), A_gmp.rowsize());
    ASSERT_EQ(QR1_gmp.colsize(), A_gmp.colsize());
    ASSERT_EQ(QR2_gmp.rowsize(), A_gmp.rowsize());
    ASSERT_EQ(QR2_gmp.colsize(), A_gmp.colsize());
    for (int i = 0; i < QR1_gmp.rowsize(); i++) {
        for (int j = 0; j < QR1_gmp.colsize(); j++) {
            ASSERT_TRUE(
                (QR1_gmp(i, j) > A_gmp(i, j) - threshold_gmp) &&
                (QR1_gmp(i, j) <
                 A_gmp(i, j) +
                     threshold_gmp));
            ASSERT_TRUE(
                (QR2_gmp(i, j) > A_gmp(i, j) - threshold_gmp) &&
                (QR2_gmp(i, j) <
                 A_gmp(i, j) +
                     threshold_gmp));
        }
    } 
}

TEST(TestQRDecompositionPivoting, TestRankReveal) {
    // declare two sample small matrices with different rank
    hdnum::DenseMatrix<double> M1({{1, 2, 3},
                                   {1, 3, -5},
                                   {1, 5, 9},
                                   {-2, 10, 15}});  // rank = 3

    hdnum::DenseMatrix<double> M2({{1, 3, 4},
                                   {1, 3, 4},
                                   {1, 1, 2},
                                   {-2, 10, 8}});  // rank = 2

    // declare two sample square matrices with different rank
    hdnum::DenseMatrix<double> M3({{1, 2, 3},
                                   {1, 3 ,-5},
                                   {2, 4, 19}});     // rank = 3

    hdnum::DenseMatrix<double> M4({{5, 2, 3},
                                   {-2, 3 ,-5},
                                   {10, 4, 6}});    // rank = 2

    // declare a sample wide matrix
    hdnum::DenseMatrix<double> M5({{1, 2, 3, 4},
                                   {-5, 9, 2, 5},
                                   {3, -9, 12, -22}});  // rank = 3 

    // declare an empty matrix
    hdnum::DenseMatrix<double> M6;    // rank = 0                       

    int rank;
    hdnum::Vector<int> p1(3);
    hdnum::Vector<int> p2(4);
    hdnum::Vector<int> p3(0);

    hdnum::qr_gram_schmidt_pivoting(M1, p1, rank);
    ASSERT_EQ(rank, 3);

    hdnum::qr_gram_schmidt_pivoting(M2, p1, rank);
    ASSERT_EQ(rank, 2);

    hdnum::qr_gram_schmidt_pivoting(M3, p1, rank);
    ASSERT_EQ(rank, 3);

    hdnum::qr_gram_schmidt_pivoting(M4, p1, rank);
    ASSERT_EQ(rank, 2);

    hdnum::qr_gram_schmidt_pivoting(M5, p2, rank);
    ASSERT_EQ(rank, 3);

    hdnum::qr_gram_schmidt_pivoting(M6, p3, rank);
    ASSERT_EQ(rank, 0);
}

TEST(TestQRDecompositionPivoting, TestSmallMatrix) {
    // declare two sample small matrices with different rank
    hdnum::DenseMatrix<double> A1({{1, 2, 3},
                                   {1, 3, -5},
                                   {1, 5, 9},
                                   {-2, 10, 15}});  // rank = 3

    hdnum::DenseMatrix<double> A2({{1, 3, 4},
                                   {1, 3, 4},
                                   {1, 1, 2},
                                   {-2, 10, 8}});  // rank = 2
    hdnum::DenseMatrix<double> Q1(A1);
    hdnum::DenseMatrix<double> Q2_temp(A2);

    int rank1;
    int rank2;
    hdnum::Vector<int> p1(3);
    hdnum::Vector<int> p2(3);
    hdnum::DenseMatrix<double> R1(hdnum::qr_gram_schmidt_pivoting(Q1, p1, rank1));
    hdnum::DenseMatrix<double> R2_temp(hdnum::qr_gram_schmidt_pivoting(Q2_temp, p2, rank2));

    // matrix A2 hasn't full rank so Q2 has the dimension mxk (k = rank(A2))
    hdnum::DenseMatrix<double> Q2(A2.rowsize(), rank2);
    for (int i = 0; i < Q2.rowsize(); i++) {
        for (int j = 0; j < Q2.colsize(); j++) {
            Q2(i, j) = Q2_temp(i, j);
        }
    }

    // matrix A2 hasn't full rank so R2 has the dimension kxn (k = rank(A2))
    hdnum::DenseMatrix<double> R2(rank2, A2.colsize());
    for (int i = 0; i < R2.rowsize(); i++) {
        for (int j = 0; j < R2.colsize(); j++) {
            R2(i, j) = R2_temp(i, j);
        }
    }

    // declare a threshold for the asserts
    double threshold = 0.00000001;

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q1.colsize(), A1.colsize());
    ASSERT_EQ(Q1.rowsize(), A1.rowsize());
    ASSERT_EQ(Q2.colsize(), rank2);
    ASSERT_EQ(Q2.rowsize(), A2.rowsize());

    // check that norm(qi) = 1 in Q1
    for (int i = 0; i < Q1.colsize(); i++) {
        double normQ1(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ1 += Q1(j, i) * Q1(j, i);
        }
        normQ1 = sqrt(normQ1);
        ASSERT_NEAR(normQ1, 1.0, threshold);
    }

    // check that norm(qi) = 1 in Q2
    for (int i = 0; i < Q2.colsize(); i++) {
        double normQ2(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ2 += Q2(j, i) * Q2(j, i);
        }
        normQ2 = sqrt(normQ2);
        ASSERT_NEAR(normQ2, 1.0, threshold);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I1(Q1.transpose() * Q1);
    hdnum::DenseMatrix<double> I2(Q2.transpose() * Q2);

    for (int i = 0; i < I1.rowsize(); i++) {
        for (int j = 0; j < I1.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I1(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I1(i, j), 0.0, threshold);
        }
    }
    for (int i = 0; i < I2.rowsize(); i++) {
        for (int j = 0; j < I2.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I2(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I2(i, j), 0.0, threshold);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R1.rowsize(), rank1);
    ASSERT_EQ(R1.colsize(), A1.colsize());
    ASSERT_EQ(R2.rowsize(), rank2);
    ASSERT_EQ(R2.colsize(), A2.colsize());

    // R is an upper triangular matrix
    for (int i = 0; i < R1.colsize(); i++) {
        for (int j = i + 1; j < R1.rowsize(); j++) {
            ASSERT_NEAR(R1(j, i), 0.0, threshold);
        }
    }
    for (int i = 0; i < R2.colsize(); i++) {
        for (int j = i + 1; j < R2.rowsize(); j++) {
            ASSERT_NEAR(R2(j, i), 0.0, threshold);
        }
    }

    // A = Q*R - maybe permutation needed
    hdnum::DenseMatrix<double> QR1(Q1 * R1);
    hdnum::DenseMatrix<double> QR2(Q2 * R2);
    hdnum::apply_permutation(QR1, p1);
    hdnum::apply_permutation(QR2, p2);
    ASSERT_EQ(QR1.rowsize(), A1.rowsize());
    ASSERT_EQ(QR1.colsize(), A1.colsize());
    ASSERT_EQ(QR2.rowsize(), A2.rowsize());
    ASSERT_EQ(QR2.colsize(), A2.colsize());
    for (int i = 0; i < QR1.rowsize(); i++) {
        for (int j = 0; j < QR1.colsize(); j++) {
            ASSERT_NEAR(QR1(i, j), A1(i, j), threshold);
        }
    }
    for (int i = 0; i < QR2.rowsize(); i++) {
        for (int j = 0; j < QR2.colsize(); j++) {
            ASSERT_NEAR(QR2(i, j), A2(i, j), threshold);
        }
    }
}

TEST(TestQRDecompositionPivoting, TestSquareMatrix) {
    // declare two sample square matrices with different rank
    hdnum::DenseMatrix<double> A1({{1, 2, 3},
                                   {1, 3 ,-5},
                                   {2, 4, 19}});     // rank = 3

    hdnum::DenseMatrix<double> A2({{5, 2, 3},
                                   {-2, 3 ,-5},
                                   {10, 4, 6}});    // rank = 2
    
    hdnum::DenseMatrix<double> Q1(A1);
    hdnum::DenseMatrix<double> Q2_temp(A2);

    int rank1;
    int rank2;
    hdnum::Vector<int> p1(3);
    hdnum::Vector<int> p2(3);
    hdnum::DenseMatrix<double> R1(hdnum::qr_gram_schmidt_pivoting(Q1, p1, rank1));
    hdnum::DenseMatrix<double> R2_temp(hdnum::qr_gram_schmidt_pivoting(Q2_temp, p2, rank2));

    // matrix A2 hasn't full rank so Q2 has the dimension mxk (k = rank(A2))
    hdnum::DenseMatrix<double> Q2(A2.rowsize(), rank2);
    for (int i = 0; i < Q2.rowsize(); i++) {
        for (int j = 0; j < Q2.colsize(); j++) {
            Q2(i, j) = Q2_temp(i, j);
        }
    }

    // matrix A2 hasn't full rank so R2 has the dimension kxn (k = rank(A2))
    hdnum::DenseMatrix<double> R2(rank2, A2.colsize());
    for (int i = 0; i < R2.rowsize(); i++) {
        for (int j = 0; j < R2.colsize(); j++) {
            R2(i, j) = R2_temp(i, j);
        }
    }

    // declare a threshold for the asserts
    double threshold = 0.00000001;

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q1.colsize(), A1.colsize());
    ASSERT_EQ(Q1.rowsize(), A1.rowsize());
    ASSERT_EQ(Q2.colsize(), rank2);
    ASSERT_EQ(Q2.rowsize(), A2.rowsize());

    // check that norm(qi) = 1 in Q1
    for (int i = 0; i < Q1.colsize(); i++) {
        double normQ1(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ1 += Q1(j, i) * Q1(j, i);
        }
        normQ1 = sqrt(normQ1);
        ASSERT_NEAR(normQ1, 1.0, threshold);
    }

    // check that norm(qi) = 1 in Q2
    for (int i = 0; i < Q2.colsize(); i++) {
        double normQ2(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ2 += Q2(j, i) * Q2(j, i);
        }
        normQ2 = sqrt(normQ2);
        ASSERT_NEAR(normQ2, 1.0, threshold);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I1(Q1.transpose() * Q1);
    hdnum::DenseMatrix<double> I2(Q2.transpose() * Q2);

    for (int i = 0; i < I1.rowsize(); i++) {
        for (int j = 0; j < I1.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I1(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I1(i, j), 0.0, threshold);
        }
    }
    for (int i = 0; i < I2.rowsize(); i++) {
        for (int j = 0; j < I2.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I2(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I2(i, j), 0.0, threshold);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R1.rowsize(), rank1);
    ASSERT_EQ(R1.colsize(), A1.colsize());
    ASSERT_EQ(R2.rowsize(), rank2);
    ASSERT_EQ(R2.colsize(), A2.colsize());

    // R is an upper triangular matrix
    for (int i = 0; i < R1.colsize(); i++) {
        for (int j = i + 1; j < R1.rowsize(); j++) {
            ASSERT_NEAR(R1(j, i), 0.0, threshold);
        }
    }
    for (int i = 0; i < R2.colsize(); i++) {
        for (int j = i + 1; j < R2.rowsize(); j++) {
            ASSERT_NEAR(R2(j, i), 0.0, threshold);
        }
    }

    // A = Q*R - maybe permutation needed
    hdnum::DenseMatrix<double> QR1(Q1 * R1);
    hdnum::DenseMatrix<double> QR2(Q2 * R2);
    hdnum::apply_permutation(QR1, p1);
    hdnum::apply_permutation(QR2, p2);
    ASSERT_EQ(QR1.rowsize(), A1.rowsize());
    ASSERT_EQ(QR1.colsize(), A1.colsize());
    ASSERT_EQ(QR2.rowsize(), A2.rowsize());
    ASSERT_EQ(QR2.colsize(), A2.colsize());
    for (int i = 0; i < QR1.rowsize(); i++) {
        for (int j = 0; j < QR1.colsize(); j++) {
            ASSERT_NEAR(QR1(i, j), A1(i, j), threshold);
        }
    }
    for (int i = 0; i < QR2.rowsize(); i++) {
        for (int j = 0; j < QR2.colsize(); j++) {
            //ASSERT_NEAR(QR2(i, j), A2(i, j), threshold);
        }
    }
}

TEST(TestQRDecompositionPivoting, TestWideMatrix) {
    // declare a sample wide matrix
    hdnum::DenseMatrix<double> A1({{1, 2, 3, 4},
                                   {-5, 9, 2, 5},
                                   {3, -9, 12, -22}});  // rank = 3 

    hdnum::DenseMatrix<double> Q1_temp(A1);

    int rank1;
    hdnum::Vector<int> p1(4);
    hdnum::DenseMatrix<double> R1_temp(hdnum::qr_gram_schmidt_pivoting(Q1_temp, p1, rank1));

    // matrix A1 hasn't full rank so Q1 has the dimension mxk (k = rank(A2))
    hdnum::DenseMatrix<double> Q1(A1.rowsize(), rank1);
    for (int i = 0; i < Q1.rowsize(); i++) {
        for (int j = 0; j < Q1.colsize(); j++) {
            Q1(i, j) = Q1_temp(i, j);
        }
    }

    // matrix A1 hasn't full rank so R1 has the dimension kxn (k = rank(A2))
    hdnum::DenseMatrix<double> R1(rank1, A1.colsize());
    for (int i = 0; i < R1.rowsize(); i++) {
        for (int j = 0; j < R1.colsize(); j++) {
            R1(i, j) = R1_temp(i, j);
        }
    }

    // declare a threshold for the asserts
    double threshold = 0.00000001;

    // check Q ∈ K^(mxn)
    ASSERT_EQ(Q1.colsize(), rank1);
    ASSERT_EQ(Q1.rowsize(), A1.rowsize());

    // check that norm(qi) = 1 in Q1
    for (int i = 0; i < Q1.colsize(); i++) {
        double normQ1(0.0);

        for (int j = 0; j < Q1.rowsize(); j++) {
            normQ1 += Q1(j, i) * Q1(j, i);
        }
        normQ1 = sqrt(normQ1);
        ASSERT_NEAR(normQ1, 1.0, threshold);
    }

    // check that Q^T * Q = I
    hdnum::DenseMatrix<double> I1(Q1.transpose() * Q1);

    for (int i = 0; i < I1.rowsize(); i++) {
        for (int j = 0; j < I1.colsize(); j++) {
            // main diagonal
            if (j == i) {
                ASSERT_NEAR(I1(i, j), 1.0, threshold);
                continue;
            }
            // other elements
            ASSERT_NEAR(I1(i, j), 0.0, threshold);
        }
    }

    // R ∈ K^(nxn)
    ASSERT_EQ(R1.rowsize(), rank1);
    ASSERT_EQ(R1.colsize(), A1.colsize());

    // R is an upper triangular matrix
    for (int i = 0; i < R1.colsize(); i++) {
        for (int j = i + 1; j < R1.rowsize(); j++) {
            ASSERT_NEAR(R1(j, i), 0.0, threshold);
        }
    }

    // A = Q*R - maybe permutation needed
    hdnum::DenseMatrix<double> QR1(Q1 * R1);
    hdnum::apply_permutation(QR1, p1);
    ASSERT_EQ(QR1.rowsize(), A1.rowsize());
    ASSERT_EQ(QR1.colsize(), A1.colsize());
    for (int i = 0; i < QR1.rowsize(); i++) {
        for (int j = 0; j < QR1.colsize(); j++) {
            ASSERT_NEAR(QR1(i, j), A1(i, j), threshold);
        }
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
