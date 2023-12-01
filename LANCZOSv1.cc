#include <iostream>
#include "hdnum.hh"

//TODO: include sparseMatrix

/**
 * UPDATE: Einbindung hdnum
 * 
 * input: sparse symmetric matrix A \in \R^(n \times n), unit 2-norm vector q_1 \in \R^n 
 * output: matrix T_k \in R^(k \times k)
 * 
 * Data structures:
 * - matrix A \in \R^(n \times n)
 * - matrix Q \in \R^(n \times k) consists of vectors q_i \in \R^n, 1 \le i \le k 
 * - scalars \alpha_i, 1 \le i \le k
 * - scalars \beta_i, 1 \le i \le k
 * - integer k, 1 \le k \le n
*/

hdnum::DenseMatrix<double> lanzcos(const hdnum::DenseMatrix<double> A, hdnum::Vector<double> q1){
    /**
     * initialisation:
     * k = 0; \beta_0 = 1, q_0 = 0, r_0 = q_1
    */
    int n = A.rowsize();
    int k = 0;
    std::vector<double> beta {1}; //beta0 = 1
    hdnum::Vector<double> q_km1(n, 0);
    hdnum::Vector<double> r(n);
    r = q1; 
    std::vector<double> alpha;

    while (beta[k] != 0 && k < n) {
        k++;
        // q_(k) = r_(k-1) / \beta_(k-1)
        auto q_k = r;
        q_k /= beta[k-1];

        // alpha_k = q_k^T * A * q_k
        hdnum::Vector<double> Aq_k(n);
        A.mv(Aq_k, q_k);
        auto alpha_k = q_k * Aq_k;
        alpha.push_back(alpha_k);
        std::cout << "alpha: " << alpha_k << std::endl;

        // r = (A - alpha_k * I) * q[k] - beta[k-1] * q[k-1]
        //(A - alpha_k * I) shortened to A(i,i) - alpha_k
        auto Amalpha = A;
        for (int i = 0; i < Amalpha.rowsize(); i++) {
            Amalpha[i][i] = Amalpha[i][i] - alpha_k;
        }
        Amalpha.mv(r, q_k);
        auto betaq = q_km1;
        betaq *= beta[k-1];
        std::cout << betaq << std::endl;
        r -= betaq; 
        std::cout << "r: " << r << std::endl;

        // beta[k] = 2-norm of r
        beta.push_back(r.two_norm());
        std::cout << "beta: " << beta[beta.size() - 1] << std::endl;

        q_km1 = q_k;
        
    }

    hdnum::DenseMatrix<double> Tk(k, k, 0);
    for(int i = 0; i < k; i++) {
        Tk[i][i] = alpha[i];
        if (i != (k-1)) {
            //Note: beta0 is not part of T_k
            Tk[i][i+1] = beta[i+1];
            Tk[i+1][i] = beta[i+1];
        }
    }

    return Tk;
}

int main(){
    //example from: 
    //https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample1_simple.cpp

    // hdnum::DenseMatrix<double> sampleMatrix(3, 3);
    // for(int i = 0; i < 3; i++) {
    //     sampleMatrix[i][i] = 2;
    // }

    hdnum::DenseMatrix<double> sampleMatrix1({{2, 1, 1},
                                              {1, 2, 1},
                                              {1, 1, 2}});
    // Its eigenvalues are {4, 1, 1}
    
    hdnum::Vector<double> sample1q1(3);
    sample1q1[0] = -0.593054;
    sample1q1[1] = 0.782547;
    sample1q1[2] = -0.189493;
    
    //std::cout << "SampleMatrix: \n" << sampleMatrix1 << std::endl;
    //std::cout << "2-norm unitvector: \n" << sample1q1 << std::endl;

    hdnum::DenseMatrix<double> Tk1 = lanzcos(sampleMatrix1, sample1q1);

    std::cout << Tk1 << std::endl;
    std::cout << "-----" << std::endl;

    hdnum::DenseMatrix<double> sampleMatrix2({{0, 1, 1},
                                              {1, 0, -1},
                                              {1, -1, 0}});
    // Its eigenvalues are {1, 1, -2}

    hdnum::Vector<double> sample2q1(3);
    sample2q1[0] = 0.818898;
    sample2q1[1] = -0.302311;
    sample2q1[2] = -0.487866;

    //std::cout << "SampleMatrix: \n" << sampleMatrix2 << std::endl;
    //std::cout << "2-norm unitvector: \n" << sample2q1 << std::endl;

    //hdnum::DenseMatrix<double> Tk2 = lanzcos(sampleMatrix2, sample2q1);

    //std::cout << Tk2 << std::endl;

    return 0;
}