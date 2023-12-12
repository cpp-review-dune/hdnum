// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_HOUSEHOLDER_HH
#define LANCZOS_HOUSEHOLDER_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>
#include "lanczos_util.hh"

/** @file
 *  @brief This file implements the Lanczos method with complete reorthogonalization using Householder matrices
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

} //namespace hdnum
#endif //LANCZOS_HOUSEHOLDER_HH