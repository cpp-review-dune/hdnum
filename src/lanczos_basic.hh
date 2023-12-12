// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_BASIC_HH
#define LANCZOS_BASIC_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>
#include "lanczos_util.hh"

/** @file
 *  @brief This file implements the basic Lanczos method
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
 * @brief basic algorithm according to book with no storage or numerical optimization
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
template<class T>
SparseMatrix<T> lanczos_basic(const SparseMatrix<T> &A, const Vector<T> &q_1){
    //init:
    int n = A.rowsize();
    int k = 0;
    std::vector<T> beta {1}; //beta_0 = 1
    Vector<T> q_km1(n, 0); //q_0 = 0
    Vector<T> r(n);
    r = q_1; //r_0 = generate_q_1
    std::vector<T> alpha;

    //TODO: additional termination condition
    while (beta[k] != 0 && k < n) {
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

    SparseMatrix<T> Tk (k,k);
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

template<class T>
SparseMatrix<T> lanczos_basic(const SparseMatrix<T> &A){
    SparseMatrix<T> Tk;
    Tk = lanczos_basic(A, generate_q_1(A.colsize())); 
    return Tk;
}

} //namespace hdnum
#endif //LANCZOS_BASIC_HH