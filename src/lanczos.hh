// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_HH
#define LANCZOS_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>

//TODO: termination condition
//TODO: eigenvalues von Tk -> ausnutzen Bandmatrix
///UPDATE: storage optimization
///UPDATE: sparseMatrix included

/** @file
 *  @brief This file implements the Lanczos method
 */

namespace hdnum {

/**
 * @brief random q_1 vector generator
 * 
 * @tparam T 
 * @param n  size of matrix \f$ A \f$
 * @return Vector<T> 2-norm unit vector \f$ \in \R^n \f$
 */
Vector<double> generate_q_1(const int &n) {
    //random T generated according to https://www.geeksforgeeks.org/generate-random-T-numbers-in-cpp/
    double lower_bound = -1;
    double upper_bound = 1;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    
    //init random vector
    Vector<double> q_1 (n);
    for (int i = 0; i < n; i++) {
        q_1[i] = unif(re);
    }

    //normalize to 2-norm
    q_1 /= q_1.two_norm();

    return q_1;
}

/**
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

/**
 * @brief computes a matrix \f$ T_k \f$ thats extremal eigenvalues are approximations 
 *        for the extreme eigenvalues of a large, sparse, symmetric matrix \f$ A \f$
 * 
 * @tparam T 
 * @param A a (large,) sparse, symmetric matrix 
 * @return SparseMatrix<T> matrix T_k
 */
template<class T>
SparseMatrix<T> lanczos(const SparseMatrix<T> &A) {
    //init
    int n = A.rowsize();
    //k = 1 equals k = 0 because our numeration starts at 0
    int k = 0;

    Vector<T> w(n);
    Vector<T> v(n);
    std::vector<T> alpha;
    std::vector<T> beta;

    //first iteration:
    //w = q_1 
    w = generate_q_1(n);
    //v = A * w
    A.mv(v, w);
    //\alpha_1 = w^T * v
    alpha.push_back(w * v);
    //v = v - \alpha_1 * w
    Vector<T> alphaw(n);
    alphaw = w;
    alphaw *= alpha[0];
    v = v - alphaw;
    //\beta_1 = ||v||_2
    beta.push_back(v.two_norm());

    // k <= 2*n for termination
    //TODO: termination
    while (beta[k] != 0.0 && k < n) {
        for (int i = 0; i < n; i++) {
            //t = w_i; w_i = v_i/\beta_k; v_i = -\beta_k * t
            T t = w[i];
            w[i] = v[i] / beta[k];
            v[i] = -1 * beta[k] * t;
        }

        //v = v + A*w
        Vector<T> Aw(n);
        A.mv(Aw, w);
        v += Aw;
        //
        k++;
        //\alpha_k = w^T * v
        alpha.push_back(w * v);
        //v = v - \alpha_k * w
        alphaw = w;
        alphaw *= alpha[0];
        v = v - alphaw;
        //\beta_k = ||v||_2
        beta.push_back(v.two_norm());
    }

    // DenseMatrix<T> Tk(k, k);
    // for(int i = 0; i < k; i++) {
    //     Tk[i][i] = alpha[i];
    //     if (i != (k-1)) {
    //         Tk[i][i+1] = beta[i];
    //         Tk[i+1][i] = beta[i];
    //     }
    // }

    SparseMatrix<T> Tk;
    auto builder = typename SparseMatrix<T>::builder(k,k);
    for (int i = 0; i < k; i++) {
        builder.addEntry(i, i, alpha[i]);
        if (i != (k-1)) {
            builder.addEntry(i, i+1, beta[i]);
            builder.addEntry(i+1, i, beta[i]);
        }
    }
    Tk = builder.build();

    return Tk;
}

/**
 * @brief basic algorithm according to book with no storage or numerical optimization
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
template<class T>
SparseMatrix<T> lanczos_basic(const SparseMatrix<T> &A){
    //init:
    int n = A.rowsize();
    int k = 0;
    std::vector<T> beta {1}; //beta_0 = 1
    Vector<T> q_km1(n, 0); //q_0 = 0
    Vector<T> r(n);
    r = generate_q_1(n); //r_0 = generate_q_1
    std::vector<T> alpha;

    //TODO: additional termination condition
    while (beta[k] != 0.0 && k < n) {
        k++;
        // q_(k) = r_(k-1) / \beta_(k-1)
        auto q_k = r;
        q_k /= beta[k-1];

        // alpha_k = q_k^T * A * q_k
        Vector<T> Aq_k(n);
        A.mv(Aq_k, q_k);
        auto alpha_k = q_k * Aq_k;
        alpha.push_back(alpha_k);

        // r = (A - alpha_k * I) * q[k] - beta[k-1] * q[k-1]
        //(A - alpha_k * I) shortened to A(i,i) - alpha_k
        auto Amalpha = A;
        for (int i = 0; i < Amalpha.rowsize(); i++) {
            Amalpha.get(i,i) = Amalpha(i,i) - alpha_k;
        }
        Amalpha.mv(r, q_k);
        auto betaq = q_km1;
        betaq *= beta[k-1];
        r -= betaq; 

        // beta[k] = 2-norm of r
        beta.push_back(r.two_norm());

        q_km1 = q_k; 
    }

    // DenseMatrix<T> Tk(k, k, 0);
    // for(int i = 0; i < k; i++) {
    //     Tk[i][i] = alpha[i];
    //     if (i != (k-1)) {
    //         //Note: beta0 is not part of T_k
    //         Tk[i][i+1] = beta[i+1];
    //         Tk[i+1][i] = beta[i+1];
    //     }
    // }

    SparseMatrix<T> Tk (k,k);
    auto builder = typename SparseMatrix<T>::builder(k,k);
    for (int i = 0; i < k; i++) {
        builder.addEntry(i, i, alpha[i]);
        if (i != (k-1)) {
            builder.addEntry(i, i+1, beta[i]);
            builder.addEntry(i+1, i, beta[i]);
        }
    }
    Tk = builder.build();

    return Tk;
}

/**
 * @brief Lanczos method with complete reorthogonalization using Householder matrices
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
//template<class T>
//SparseMatrix<T> lanczos_householder(const SparseMatrix<T> &A) {
    /*
     * Algorithm:
     * 
     * r_0 = q_1 (given unit vector)
     * Determine Householder H_0 s.t. H_0 * r_0 = e_1
     * for (k = 1 .. n-1) {
     *      a_k = q_k^T * A * q_k
     *      r_k = (A-alpha_k * I) * q_k - beta_(k-1) * q_(k-1)
     *      w = (H_(k-1) * ... * H_0) * r_k
     *      Determine Householder H_k s.t. H_k * w = [w_1, ..., w_k, beta_k, 0,...,0]^T
     *      q_(k+1) = H_0 * .... * H_k * e_(k+1)
     * }
     */

//     SparseMatrix<T> Tk(k,k);
//     return T_k;
// }


/**
 * @brief OTHER ALGORITHM ATTEMPT
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
template<class T>
SparseMatrix<T> lanczos_complete(const SparseMatrix<T> &A) {
    /* Algorithm:
     * q = x / ||x||
     * Q_1 = [q]
     * r = A * q
     * alpha_1 = q * r
     * r = r - alpha_1 * q
     * beta_1 = ||r||
     * j = 2
     * 
     * while(beta_j != 0.0) {
     *      v = q
     *      q = r/beta_(j-1)
     *      Q_j = [Q_(j-1), q]
     *      r = A * q - beta_(j-1) * v
     *      alpha_j = q * r
     *      r = r - alpha_j * q
     *      r = r - Q * (Q * r)
     *      beta_j = ||r||
     *      j++
     * } 
     */

    //init:
    int n = A.rowsize();
    //k = 2
    int k = 1;

    std::vector<T> alpha;
    std::vector<T> beta;
    Vector<T> q(n);
    DenseMatrix<T> Q;
    Vector<T> r(n);
    Vector<T> v(n);

    //q_1 = x / ||x||; Q = [q]
    q = generate_q_1(n);
    Q.transpose(); Q.addNewRow(q); Q.transpose();

    //first iteration
    //r = A * q
    A.mv(r, q);
    //alpha_1 = q * r
    alpha.push_back(q * r);
    //r = r - alpha_1 * q
    r.update(-1 * alpha.back(), q);
    //beta_1 = ||r||
    beta.push_back(r.two_norm());
    //std::cout << r << std::endl;
    
    while(beta[k] != 0.0 && k < n) {
        //v = q_(k-1)
        auto q_km1 = q;

        //q = r/beta_(j-1)
        q = r;
        q /= beta[k-1];

        //Q_j = [Q_(j-1), q]
        Q.transpose(); Q.addNewRow(q); Q.transpose();

        //r = A * q - beta_(j-1) * v
        A.mv(r, q);
        r.update(-1 * beta.back(), q_km1);

        //alpha_j = q * r
        alpha.push_back(q * r);

        //r = r - alpha_j * q
        r.update(-1 * alpha.back(), q);

        //reorthogonalization
        //r = r - Q * (Q^T * r)
        Vector<T> QQr(n);
        (Q.transpose()).mv(QQr, r);
        Q.mv(QQr, QQr);
        r -= QQr;
        //std::cout << r << std::endl;

        //beta_j = ||r||
        beta.push_back(r.two_norm());

        k++;
    }

    SparseMatrix<T> Tk;
    auto builder = typename SparseMatrix<T>::builder(k,k);
    for (int i = 0; i < k; i++) {
        builder.addEntry(i, i, alpha[i]);
        if (i != (k-1)) {
            builder.addEntry(i, i+1, beta[i]);
            builder.addEntry(i+1, i, beta[i]);
        }
    }
    Tk = builder.build();

    return Tk;
}

/**
 * @brief 
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
// template<class T>
// SparseMatrix<T> lanczos_selective(const SparseMatrix<T> &A){
//     SparseMatrix<T> Tk(k,k);
//     return Tk;
// }

} //namespace hdnum
#endif //LANCZOS_HH