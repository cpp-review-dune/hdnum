// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_SRO_HH
#define LANCZOS_SRO_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>
#include "lanczos_util.hh"

/** @file
 *  @brief This file implements the Lanczos method
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
 * @brief 
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
// template<class T>
// SparseMatrix<T> lanczos_sro(const SparseMatrix<T> &A){
//     SparseMatrix<T> Tk(k,k);
//     return Tk;
// }

} //namespace hdnum
#endif //LANCZOS_SRO_HH