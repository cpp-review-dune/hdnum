// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_CRO_HH
#define LANCZOS_CRO_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>
#include "lanczos_util.hh"

/** @file
 *  @brief This file implements the Lanczos method with complete reorthogonalization
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
 * @brief OTHER ALGORITHM ATTEMPT
 * 
 * @tparam T 
 * @param A 
 * @return SparseMatrix<T> 
 */
template<class T>
SparseMatrix<T> lanczos_cro(const SparseMatrix<T> &A) {
    //init:
    int n = A.rowsize();
    int k = 0;

    std::vector<T> alpha;
    std::vector<T> beta;
    Vector<T> q(n);
    DenseMatrix<T> Q;
    DenseMatrix<T> QT;
    Vector<T> r(n);
    Vector<T> v(n);

    //q = x/||x||; Q = [q]
    q = generate_q_1(n);
    QT = Q.transpose(); QT.addNewRow(q); Q = QT.transpose();

    //first iteration
    //r = A * q
    A.mv(r, q);
    //alpha_1 = q * r
    alpha.push_back(q * r);
    //r = r - alpha_1 * q
    r.update(-1 * alpha.back(), q);
    //beta_1 = ||r||
    beta.push_back(r.two_norm());
    
    while(std::abs(beta[k]) > 1e-9 && k < n) {
        k++;
        //v = q_(k-1)
        auto q_km1 = q;

        //q = r/beta_(j-1)
        q = r;
        q /= beta[k-1];

        //Q_j = [Q_(j-1), q]
        QT = Q.transpose(); QT.addNewRow(q); Q = QT.transpose();

        //r = A * q - beta_(j-1) * v
        A.mv(r, q);
        r.update(-1 * beta.back(), q_km1);

        //alpha_j = q * r
        alpha.push_back(q * r);

        //r = r - alpha_j * q
        r.update(-1 * alpha.back(), q);

        //reorthogonalization
        //r = r - Q * (Q^T * r)
        Vector<T> QTr (Q.colsize());
        QT.mv(QTr, r);
        Vector<T> QQTr (n);
        Q.mv(QQTr, QTr);
        r -= QQTr;

        //beta_j = ||r||
        beta.push_back(r.two_norm());

        //test:
        DenseMatrix<T> Tk_d(k, k);
        for(int i = 0; i < k; i++) {
            Tk_d[i][i] = alpha[i];
            if (i != (k-1)) {
                Tk_d[i][i+1] = beta[i];
                Tk_d[i+1][i] = beta[i];
            }
        }
        std::vector<T> real;
        std::vector<T> imag;

        eigenvalues_qr_algorithm_givens(Tk_d, real, imag);
        for (auto i : real) {
            if (i == 200.000) {
                std::cout << "k:" << k << std::endl;
            }
        }
        std::cout << beta[k] << std::endl;
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

} //namespace hdnum
#endif //LANCZOS_CRO_HH