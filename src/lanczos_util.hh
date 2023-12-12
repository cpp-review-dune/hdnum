// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_UTIL_HH
#define LANCZOS_UTIL_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>

// #include "lanczos_basic.hh"
// #include "lanczos_optimized.hh"
// #include "lanczos_householder.hh"
// #include "lanczos_cro.hh"
// #include "lanczos_sro.hh"

//TODO: termination condition
//TODO: eigenvalues von Tk -> ausnutzen Bandmatrix
//TODO: test application without randomized q_1

///UPDATE: storage optimization
///UPDATE: sparseMatrix included
///UPDATE: first CRO method implemented

/** @file
 *  @brief This file implements the Lanczos method
 */

namespace hdnum {

/**
 * @brief random q_1 vector generator
 * 
 * @param n  size of matrix \f$ A \f$
 * @return Vector<T> 2-norm unit vector \f$ \in \R^n \f$
 */
Vector<double> generate_q_1(const size_t &n) {
    //random double generated according to https://www.geeksforgeeks.org/generate-random-T-numbers-in-cpp/
    double lower_bound = -1;
    double upper_bound = 1;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    
    //init random vector
    Vector<double> q_1 (n);
    for (size_t i = 0; i < n; i++) {
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

} //namespace hdnum
#endif //LANCZOS_UTIL_HH