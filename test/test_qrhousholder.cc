#include <ctime>

#include "../hdnum.hh"
#include "gtest/gtest.h"

// test R and Vi vectors for quadratic matrix

TEST(TestQRhousholder, TestPostconditionsqrhousholderRandVi) {
    // sample matrix where we know the expected results
    hdnum::DenseMatrix<double> A = {
        {1, -1, 4}, {1, 4, -2}, {1, 4, 2}, {1, -1, 0}};
    hdnum::Vector<double> v(A.colsize(), 0);
    double threshold = 0.000001;
    // test if R and V_{i} correct after calling the function
    hdnum::qrhousholder(A, v);
    ASSERT_NEAR(A(0, 0), 3.0, threshold);
    ASSERT_NEAR(A(1, 0), 1.0, threshold);
    ASSERT_NEAR(A(2, 0), 1.0, threshold);
    ASSERT_NEAR(A(3, 0), 1.0, threshold);
    ASSERT_NEAR(A(0, 1), -3.0, threshold);
    ASSERT_NEAR(A(1, 1), 8.333333, threshold);
    ASSERT_NEAR(A(2, 1), 3.333333, threshold);
    ASSERT_NEAR(A(3, 1), -1.666667, threshold);
    ASSERT_NEAR(A(0, 2), -2.0, threshold);
    ASSERT_NEAR(A(1, 2), 2.0, threshold);
    ASSERT_NEAR(A(2, 2), 6.4, threshold);
    ASSERT_NEAR(A(3, 2), -3.2, threshold);
    ASSERT_NEAR(v[0], -2.0, threshold);
    ASSERT_NEAR(v[1], -5.0, threshold);
    ASSERT_NEAR(v[2], -4.0, threshold);
}

// test R and Vi vectors for quadratic matrix

TEST(TestQRhousholder,
     TestQRhousholder_Test_Postconditions_qrhousholder_quadratisch) {
    // sample matrix where we know the expected results
    hdnum::DenseMatrix<double> A = {{1, -1, 4}, {1, 4, -2}, {1, 4, 2}};
    hdnum::Vector<double> V(A.colsize(), 0);
    double threshold = 0.0001;
    // test if R and V_{i} correct after calling the function
    qrhousholder(A, V);
    ASSERT_NEAR(A(0, 0), 2.7321, threshold);
    ASSERT_NEAR(A(1, 0), 1.0, threshold);
    ASSERT_NEAR(A(2, 0), 1.0, threshold);
    ASSERT_NEAR(A(0, 1), -4.0415, threshold);
    ASSERT_NEAR(A(1, 1), 6.9692, threshold);
    ASSERT_NEAR(A(2, 1), 2.8868, threshold);
    ASSERT_NEAR(A(0, 2), -2.3094, threshold);
    ASSERT_NEAR(A(1, 2), 3.266, threshold);
    ASSERT_NEAR(A(2, 2), 5.6569, threshold);
    ASSERT_NEAR(V[0], -1.7321, threshold);
    ASSERT_NEAR(V[1], -4.0825, threshold);
    ASSERT_NEAR(V[2], -2.8284, threshold);
}

// test Q = Qt for a random matrix
TEST(TestQRhousholder, TestPostconditionsqrhousholderexplizitQ) {
    // create a random size small matrix
    int m;
    int n;
    double threshold = 0.0000001;
    do {
        srand(time(NULL));
        m = (rand() % 18) + 2;  // [2, ... ,20]
        n = (rand() % 18) + 2;  // [2, ... ,20]
    } while (m <= n);
    hdnum::DenseMatrix<double> A(m, n);
    hdnum::Vector<double> v(A.colsize(), 0);
    // fill it with random elements
    for (int i = 0; i < A.rowsize(); i++) {
        for (int j = 0; j < A.colsize(); j++) {
            int x = (rand() % 200) - 100;  // [-100, ... ,100]
            A(i, j) = x;
        }
    }
    hdnum::DenseMatrix<double> I(m, m);
    hdnum::identity(I);
    hdnum::DenseMatrix<double> q = hdnum::qrhousholderexplizitQ(A, v, false);
    hdnum::DenseMatrix<double> res = q * q.transpose();
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            if (j == i) {
                ASSERT_NEAR(res(i, j), 1.0, threshold);
            }
            if (j != i) {
                ASSERT_NEAR(res(i, j), 0.0, threshold);
            }
        }
    }
}

// test Q * R = A for m > n
TEST(TestQRhousholder, TestPostconditionsqrhousholderexplizitQ_results) {
    // create a random size small matrix
    int m;
    int n;
    double threshold = 0.0000001;
    do {
        srand(time(NULL));
        m = (rand() % 4) + 2;  // [2, ... ,20]
        n = (rand() % 4) + 2;  // [2, ... ,20]
    } while (m <= n);
    hdnum::DenseMatrix<double> A(m, n);
    hdnum::Vector<double> v(A.colsize(), 0);
    // fill it with random elements
    for (int i = 0; i < A.rowsize(); i++) {
        for (int j = 0; j < A.colsize(); j++) {
            int x = (rand() % 10) - 5;  // [-100, ... ,100]
            A(i, j) = x;
        }
    }
    hdnum::DenseMatrix<double> orginal(A);
    hdnum::DenseMatrix<double> I(m, m);
    hdnum::identity(I);
    hdnum::DenseMatrix<double> q = hdnum::qrhousholderexplizitQ(A, v, false);

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                A(i, j) = v[i];
                continue;
            }
            if (i > j) {
                A(i, j) = 0;
            }
        }
    }
    hdnum::DenseMatrix<double> res(q * A);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            ASSERT_NEAR(res(i, j), orginal(i, j), threshold);
        }
    }
}
TEST(TestQRhousholder,
     TestPostconditionsqrhousholderexplizitQ_results_quadrat) {
    // create a random size small matrix
    int m;
    int n;
    double threshold = 0.0000001;
    do {
        srand(time(NULL));
        m = (rand() % 4) + 2;  // [2, ... ,20]
        n = (rand() % 4) + 2;  // [2, ... ,20]
    } while (m != n);
    hdnum::DenseMatrix<double> A(m, n);
    hdnum::Vector<double> v(A.colsize(), 0);
    // fill it with random elements
    for (int i = 0; i < A.rowsize(); i++) {
        for (int j = 0; j < A.colsize(); j++) {
            int x = (rand() % 10) - 5;  // [-100, ... ,100]
            A(i, j) = x;
        }
    }
    hdnum::DenseMatrix<double> orginal(A);
    hdnum::DenseMatrix<double> I(m, m);
    hdnum::identity(I);
    hdnum::DenseMatrix<double> q = hdnum::qrhousholderexplizitQ(A, v, false);

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                A(i, j) = v[i];
                continue;
            }
            if (i > j) {
                A(i, j) = 0;
            }
        }
    }
    hdnum::DenseMatrix<double> res(q * A);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            ASSERT_NEAR(res(i, j), orginal(i, j), threshold);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}