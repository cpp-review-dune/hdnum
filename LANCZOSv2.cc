#include <iostream>
#include "hdnum.hh"

//TODO: include sparseMatrix

/**
 * UPDATE: storage optimization
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
 * 
 * randomize q_1 vector
 * 
*/

hdnum::DenseMatrix<double> lanzcos(const hdnum::DenseMatrix<double> A, hdnum::Vector<double> q1){
    /**
     * new initialisation:
     * The first iteration of the algorithm is done in the initialisation.
     * w = q_1; v = A * w; \alpha_1 = w^T * v; v = v - \alpha_1 * w; \beta_1 = ||v||_2; k = 1
    */
    int n = A.rowsize();
    //k = 1 equals k = 0 because our numeration starts with 0
    int k = 0;

    //init
    hdnum::Vector<double> w(n);
    hdnum::Vector<double> v(n);
    std::vector<double> alpha;
    std::vector<double> beta;

    //first iteration:
    //w = q_1; v = A * w; \alpha_1 = w^T * v; v = v - \alpha_1 * w; \beta_1 = ||v||_2
    w = q1;
    //
    A.mv(v, w);
    //
    alpha.push_back(w * v);
    std::cout << "alpha: " << alpha[alpha.size() - 1] << std::endl;
    //
    hdnum::Vector<double> alphaw(n);
    alphaw = w;
    alphaw *= alpha[0];
    std::cout << alphaw << std::endl;
    v = v - alphaw;
    std::cout << "r: " << v << std::endl;
    //
    beta.push_back(v.two_norm());
    std::cout << "beta: " << beta[beta.size() - 1] << std::endl;

    // k <= 2*n for termination
    while (beta[k] != 0 && k < (2*n)) {
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
        std::cout << "alpha: " << alpha[alpha.size()-1] << std::endl;

        //v = v - \alpha_k * w
        alphaw = w;
        alphaw *= alpha[0];
        std::cout << alphaw << std::endl;
        v = v - alphaw;
        std::cout << "r: " << v << std::endl;

        //\beta_k = ||v||_2
        beta.push_back(v.two_norm());
        std::cout << "beta: " << beta[beta.size() - 1] << std::endl;
    }

    hdnum::DenseMatrix<double> Tk(k, k, 0);
    for(int i = 0; i < k; i++) {
        Tk[i][i] = alpha[i];
        if (i != (k-1)) {
            Tk[i][i+1] = beta[i];
            Tk[i+1][i] = beta[i];
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