// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/*
 * File:   qrhousholder
 * Author: Ahmad Fadl <abohmaid@windowslive.com>
 *
 * Created on August 25, 2020
 */

#ifndef HDNUM_QRHOUSHOLDER_HH
#define HDNUM_QRHOUSHOLDER_HH
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "densematrix.hh"
#include "vector.hh"
/** @file
 *  @brief This file implements QR decomposition using housholder transformation
 */
namespace hdnum {
template <class REAL>
DenseMatrix<REAL> creat_I_matrix(size_t n) {
    DenseMatrix<REAL> res(n, n, 0);
    for (size_t i = 0; i < n; i++) {
        res(i, i) = 1;
    }
    return res;
}

/*! \brief Function that return the sign of a number
 */
template <typename REAL>
size_t sgn(REAL val) {
    return (REAL(0) < val) - (val < REAL(0));
}
/*! @brief Funktion that calculate the QR decoposition in place
the elements of A will be replaced with the elements of v_{i}vectors and the
upper diagonals elements of R and the diagonal elements of R will be saved in
vectro v.


\tparam A the Matrix
\tparam v oa vector of hdnum::Vector
*/
template <class REAL>
void qrhousholder(DenseMatrix<REAL>& A, hdnum::Vector<REAL>& v) {
    auto m = A.rowsize();
    auto n = A.colsize();
    for (size_t j = 0; j < n; j++) {
        REAL s = 0;
        for (size_t i = j; i < m; i++) {
            s = s + pow(A(i, j), 2);
        }
        s = sqrt(s);
        v[j] = (-1.0) * sgn(A(j, j)) * s;
        REAL fak = sqrt(s * (s + std::abs(A(j, j))));
        A(j, j) = A(j, j) - v[j];
        for (size_t k = j; k < m; k++) {
            A(k, j) = A(k, j) / fak;
        }
        for (size_t i = j + 1; i < n; i++) {
            s = 0;
            for (size_t k = j; k < m; k++) {
                s = s + (A(k, j) * A(k, i));
            }
            for (size_t k = j; k < m; k++) {
                A(k, i) = A(k, i) - (A(k, j) * s);
            }
        }
        // normalize the vi vectors again
        for (size_t i = m; i >= 0; i--) {
            A(i, j) = A(i, j) * fak;
            if (i == j) {
                break;
            }
        }
    }
}
/*! @brief Funktion that calculate the QR decoposition in place and return Q
the elements of A will be replaced with the elements of v_{i}vectors and the
upper diagonals elements of R and the diagonal elements of R will be saved in
vectro v.


\tparam A the Matrix
\tparam v oa vector of hdnum::Vector
\return Q matrix
*/
template <class REAL>
DenseMatrix<REAL> qrhousholderexplizitQ(DenseMatrix<REAL>& A,
                                        hdnum::Vector<REAL>& v,
                                        bool show_Hi = false) {
    auto m = A.rowsize();
    auto n = A.colsize();
    auto I = creat_I_matrix<REAL>(std::max(m, n));

    DenseMatrix<REAL> Q(m, m, 0);
    for (size_t j = 0; j < n; j++) {
        REAL s = 0;
        for (size_t i = j; i < m; i++) {
            s = s + pow(A(i, j), 2);
        }
        s = sqrt(s);
        v[j] = (-1.0) * sgn(A(j, j)) * s;
        REAL fak = sqrt(s * (s + std::abs(A(j, j))));
        A(j, j) = A(j, j) - v[j];
        for (size_t k = j; k < m; k++) {
            A(k, j) = A(k, j) / fak;
        }
        for (size_t i = j + 1; i < n; i++) {
            s = 0;
            for (size_t k = j; k < m; k++) {
                s = s + (A(k, j) * A(k, i));
            }
            for (size_t k = j; k < m; k++) {
                A(k, i) = A(k, i) - (A(k, j) * s);
            }
        }
        // normalize the vi vectors again
        for (size_t i = m; i >= 0; i--) {
            A(i, j) = A(i, j) * fak;
            if (i == j) {
                break;
            }
        }
    }
    // create qi and multiply them
    if (m >= n) {
        for (size_t j = 0; j < n; j++) {
            DenseMatrix<REAL> TempQ(m, m, 0.0);
            DenseMatrix<REAL> v1(m, 1, 0.0);
            DenseMatrix<REAL> v1t(1, m, 0.0);
            hdnum::Vector<double> v__i(m, 0);
            for (size_t i = 0; i < m; i++) {
                if (i < j) {
                    v1(i, 0) = 0;

                    v__i[i] = 0;
                    continue;
                }
                v1(i, 0) = A(i, j);

                v__i[i] = A(i, j);
            }
            v1t = v1.transpose();

            TempQ = (v1 * v1t);

            TempQ *= (-2.0);

            TempQ /= v__i.two_norm_2();

            TempQ += I;
            if (show_Hi) {
                std::cout << "H[" << j + 1 << "]" << TempQ;
            }
            if (j == 0) {
                Q = TempQ;
            }
            if (j > 0) {
                Q = Q * TempQ;
            }
        }
    }
    if (n > m) {
        for (size_t j = 0; j < m; j++) {
            DenseMatrix<REAL> TempQ(m, m, 0.0);
            DenseMatrix<REAL> v1(m, 1, 0.0);
            DenseMatrix<REAL> v1t(1, m, 0.0);
            hdnum::Vector<double> v__i(m, 0);
            for (size_t i = 0; i < m; i++) {
                if (i < j) {
                    v1(i, 0) = 0;

                    v__i[i] = 0;
                    continue;
                }
                v1(i, 0) = A(i, j);

                v__i[i] = A(i, j);
            }
            v1t = v1.transpose();

            TempQ = (v1 * v1t);

            TempQ *= (-2.0);

            TempQ /= v__i.two_norm_2();

            TempQ += I;
            if (show_Hi) {
                std::cout << "H[" << j + 1 << "]" << TempQ;
            }
            if (j == 0) {
                Q = TempQ;
            }
            if (j > 0) {
                Q = Q * TempQ;
            }
        }
    }
    return Q;
}
}  // namespace hdnum
#endif