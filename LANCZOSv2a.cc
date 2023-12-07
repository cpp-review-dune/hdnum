// -*- tab-width: 4; indent-tabs-mode: nil -*-
#include <iostream>
#include "hdnum.hh"
#include <random>

//TODO: termination condition

/**
 * @brief random q_1 vector construction
 * 
 * @param n size of matrix \f$ A \f$
 * @return 2-norm unit vector \f$ \in \R^n \f$
 */
hdnum::Vector<double> q_1(const double n) {
    //random double generated according to https://www.geeksforgeeks.org/generate-random-double-numbers-in-cpp/
    double lower_bound = -1;
    double upper_bound = 1;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    
    //init random vector
    hdnum::Vector<double> q_1 (n);
    for (int i = 0; i < n; i++) {
        q_1[i] = unif(re);
    }
    
    std::cout << q_1 << std::endl;

    //normalize to 2-norm
    q_1 /= q_1.two_norm();

    std::cout << q_1 << std::endl;

    return q_1;
}

/**
 * UPDATE: storage optimization
 * UPDATE: sparseMatrix included
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

hdnum::SparseMatrix<double> lanzcos(const hdnum::SparseMatrix<double> &A) {
    /**
     * new initialisation:
     * The first iteration of the algorithm is done in the initialisation.
     * w = q_1; v = A * w; \alpha_1 = w^T * v; v = v - \alpha_1 * w; \beta_1 = ||v||_2; k = 1
    */
    int const n = A.rowsize();
    //k = 1 equals k = 0 because our numeration starts at 0
    int k = 0;

    //init
    hdnum::Vector<double> w(n);
    hdnum::Vector<double> v(n);
    std::vector<double> alpha;
    std::vector<double> beta;
   
    //first iteration:
    //w = q_1; v = A * w; \alpha_1 = w^T * v; v = v - \alpha_1 * w; \beta_1 = ||v||_2
    w = q_1(n);
    //
    A.mv(v, w);
    //
    alpha.push_back(w * v);
    //
    hdnum::Vector<double> alphaw(n);
    alphaw = w;
    alphaw *= alpha[0];
    //std::cout << alphaw << std::endl;
    v = v - alphaw;
    //std::cout << "r: " << v << std::endl;
    //
    beta.push_back(v.two_norm());
    //std::cout << "beta: " << beta[beta.size() - 1] << std::endl;

    // k <= 2*n for termination
    // && k < (n*n)
    while (beta[k] != 0  && k < n) {
        for (int i = 0; i < n; i++) {
            //t = w_i; w_i = v_i/\beta_k; v_i = - \beta_k * t
            double t = w[i];
            w[i] = v[i] / beta[k];
            v[i] = -1 * beta[k] * t;
        }

        //v = v + A*w
        hdnum::Vector<double> Aw(n);
        A.mv(Aw, w);
        v += Aw;

        k++;

        //\alpha_k = w^T * v
        alpha.push_back(w * v);
        //std::cout << "alpha: " << alpha[alpha.size()-1] << std::endl;

        //v = v - \alpha_k * w
        alphaw = w;
        alphaw *= alpha[0];
        //std::cout << alphaw << std::endl;
        v = v - alphaw;
        //std::cout << "r: " << v << std::endl;

        //\beta_k = ||v||_2
        beta.push_back(v.two_norm());
        //std::cout << "beta: " << beta[beta.size() - 1] << std::endl;
    }

    // hdnum::DenseMatrix<double> Tk(k, k);
    // for(int i = 0; i < k; i++) {
    //     Tk[i][i] = alpha[i];
    //     if (i != (k-1)) {
    //         Tk[i][i+1] = beta[i];
    //         Tk[i+1][i] = beta[i];
    //     }
    // }

    hdnum::SparseMatrix<double> Tk(k,k);
    auto builder = hdnum::SparseMatrix<double>::builder(k, k);
    for (int i = 0; i < k; i++) {
        builder.addEntry(i, i, alpha[i]);
        if (i != (k-1)) {
            builder.addEntry(i, i+1, beta[i]);
            builder.addEntry(i+1, i, beta[i]);
        }
    }
    Tk = builder.build();

    std::cout << Tk << std::endl;

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

    //hdnum::DenseMatrix<double> Tk1 = lanzcos(sampleMatrix1, sample1q1);

    //std::cout << Tk1 << std::endl;
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

    hdnum::SparseMatrix<double> test {};
    hdnum::readMatrixFromFile("test/matrix_market_files/lanczos_example_bcsstm01.mtx", test);

    //std::cout << test << std::endl;

    hdnum::SparseMatrix<double> Tk2 = lanzcos(test);
    std::cout << Tk2 << std::endl;

    return 0;
}