#include <ctime>

#include "../hdnum.hh"
#include "gtest/gtest.h"
TEST(TestQRhousholder, TestPostconditionsqrhousholderRandV_{i}){
    //sample matrix where we know the expected results
    hdnum::DenseMatrix<double> A={{1, -1, 4},{1 ,4, -2},{1 ,4, 2},{1 ,-1 ,0}};
    std::vector<double>v(A.colsize(),0) ;
    double threshold = 0.00000001;
    //test if R and V_{i} correct after calling the function
    qrhousholder(A,v);
    ASSERT_NEAR(A(0,0), 3.0, threshold);
    ASSERT_NEAR(A(1,0), 1.0, threshold);
    ASSERT_NEAR(A(2,0), 1.0, threshold);
    ASSERT_NEAR(A(3,0), 1.0, threshold);
    ASSERT_NEAR(A(0,1), -3.0, threshold);
    ASSERT_NEAR(A(1,1), 8.0, threshold);
    ASSERT_NEAR(A(2,1), 3.333333, threshold);
    ASSERT_NEAR(A(3,1), -1.666667, threshold);
    ASSERT_NEAR(A(0,2), -2.0, threshold);
    ASSERT_NEAR(A(1,2), 2.0, threshold);
    ASSERT_NEAR(A(2,2), 6.4, threshold);
    ASSERT_NEAR(A(3,2), -3.2, threshold);
    ASSERT_NEAR(V[0], -2.0, threshold);
    ASSERT_NEAR(V[1], -5.0, threshold);
    ASSERT_NEAR(V[2], -4.0, threshold);
}
TEST(TestQRhousholder, TestPostconditionsqrhousholderwithq){
// create a random size small matrix
    int m;
    int n;
    do {
        srand(time(NULL));
        m = (rand() % 18) + 2;  // [2, ... ,20]
        n = (rand() % 18) + 2;  // [2, ... ,20]
    } while (m <= n);
    hdnum::DenseMatrix<double> Q(m, n);

    // fill it with random elements
    for (int i = 0; i < Q.rowsize(); i++) {
        for (int j = 0; j < Q.colsize(); j++) {
            int x = (rand() % 200) - 100;  // [-100, ... ,100]
            Q(i, j) = x;
        }
    }
    auto q = hdnum::qrhousholderexplizitQ(A, v,false);
    ASSERT_TRUE((q * q.transpose())==creat_I_matrix(m));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}