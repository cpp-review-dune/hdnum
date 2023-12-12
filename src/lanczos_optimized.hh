// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_OPTIMIZED_HH
#define LANCZOS_OPTIMIZED_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>
#include "lanczos_util.hh"

/** @file
 *  @brief This file implements the storage optimized Lanczos method
 */

namespace hdnum {

// Vector<double> generate_q_1(const size_t &n) {
//     //random double generated according to https://www.geeksforgeeks.org/generate-random-T-numbers-in-cpp/
//     double lower_bound = -1;
//     double upper_bound = 1;
//     std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
//     std::default_random_engine re;
    
//     //init random vector
//     Vector<double> q_1 (n);
//     for (size_t i = 0; i < n; i++) {
//         q_1[i] = unif(re);
//     }

//     //normalize to 2-norm
//     q_1 /= q_1.two_norm();

//     return q_1;
// }

/**
 * @brief computes a matrix \f$ T_k \f$ thats extremal eigenvalues are approximations 
 *        for the extreme eigenvalues of a large, sparse, symmetric matrix \f$ A \f$
 * 
 * @tparam T 
 * @param A a (large,) sparse, symmetric matrix 
 * @return SparseMatrix<T> matrix T_k
 */
template<class T>
SparseMatrix<T> lanczos_optimized(const SparseMatrix<T> &A) {
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
    while (beta[k] != 0 && k < n) {
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
    for (size_t i = 0; i < k; i++) {
        builder.addEntry(i, i, alpha[i]);
        if (i != (k-1)) {
            builder.addEntry(i, i+1, beta[i]);
            builder.addEntry(i+1, i, beta[i]);
        }
    }
    Tk = builder.build();
    return Tk;
}

} //namespace hdnum
#endif //LANCZOS_OPTIMIZED_HH