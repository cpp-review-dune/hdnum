// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_QR_HH
#define HDNUM_QR_HH

#include <cmath>
#include <utility>

#include "densematrix.hh"
#include "vector.hh"

/** @file
 *  @brief This file implements QR decomposition using Gram-Schmidt method
 */

namespace hdnum {
//! computes orthonormal basis of Im(A) using classical Gram-Schmidt
template <class T>
DenseMatrix<T> gram_schmidt(const DenseMatrix<T>& A) {
    DenseMatrix<T> Q(A);

    // for all columns except the first
    for (int k = 1; k < Q.colsize(); k++) {
        // orthogonalize column k against all previous
        for (int j = 0; j < k; j++) {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i = 0; i < Q.rowsize(); i++) {
                sum_nom += A[i][k] * Q[i][j];
                sum_denom += Q[i][j] * Q[i][j];
            }
            // modify
            T alpha = sum_nom / sum_denom;
            for (int i = 0; i < Q.rowsize(); i++) Q[i][k] -= alpha * Q[i][j];
        }
    }
    for (int j = 0; j < Q.colsize(); j++) {
        // compute norm of column j
        T sum(0.0);
        for (int i = 0; i < Q.rowsize(); i++) sum += Q[i][j] * Q[i][j];
        sum = sqrt(sum);
        // scale
        for (int i = 0; i < Q.rowsize(); i++) Q[i][j] = Q[i][j] / sum;
    }
    return Q;
}

//! computes orthonormal basis of Im(A) using modified Gram-Schmidt
template <class T>
DenseMatrix<T> modified_gram_schmidt(const DenseMatrix<T>& A) {
    DenseMatrix<T> Q(A);

    for (int k = 0; k < Q.colsize(); k++) {
        // modify all later columns with column k
        for (int j = k + 1; j < Q.colsize(); j++) {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i = 0; i < Q.rowsize(); i++) {
                sum_nom += Q[i][j] * Q[i][k];
                sum_denom += Q[i][k] * Q[i][k];
            }
            // modify
            T alpha = sum_nom / sum_denom;
            for (int i = 0; i < Q.rowsize(); i++) Q[i][j] -= alpha * Q[i][k];
        }
    }
    for (int j = 0; j < Q.colsize(); j++) {
        // compute norm of column j
        T sum(0.0);
        for (int i = 0; i < Q.rowsize(); i++) sum += Q[i][j] * Q[i][j];
        sum = sqrt(sum);
        // scale
        for (int i = 0; i < Q.rowsize(); i++) Q[i][j] = Q[i][j] / sum;
    }
    return Q;
}

//! computes qr decomposition using modified Gram-Schmidt - works only with small (m>n) and square matrices
template <class T>
DenseMatrix<T> qr_gram_schmidt_simple(DenseMatrix<T>& Q) {
    // save matrix A, before it's replaced with Q
    DenseMatrix<T> A(Q);

    // create matrix R
    DenseMatrix<T> R(Q.colsize(), Q.colsize());

    // start orthogonalizing
    for (int k = 0; k < Q.colsize(); k++) {
        // modify all later columns with column k
        for (int j = k + 1; j < Q.colsize(); j++) {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i = 0; i < Q.rowsize(); i++) {
                sum_nom += Q(i, j) * Q(i, k);
                sum_denom += Q(i, k) * Q(i, k);
            }

            T alpha = sum_nom / sum_denom;
            for (int i = 0; i < Q.rowsize(); i++) Q(i, j) -= alpha * Q(i, k);
        }
    }

    // add values to R, except main diagonal
    for (int i = 1; i < R.colsize(); i++) {
        for (int j = 0; j < i; j++) {
            T sum_nom(0.0);
            T sum_l2nom(0.0);
            for (int k = 0; k < Q.rowsize(); k++) {
                sum_nom += A(k, i) * Q(k, j);
                sum_l2nom += Q(k, j) * Q(k, j);
            }
            sum_l2nom = sqrt(sum_l2nom);
            // add element
            R(j, i) = sum_nom / sum_l2nom;
        }
    }

    // add missing values and scale
    for (int j = 0; j < Q.colsize(); j++) {
        // compute norm of column j
        T sum(0.0);
        for (int i = 0; i < Q.rowsize(); i++) sum += Q(i, j) * Q(i, j);
        sum = sqrt(sum);
        // add main diagonal to R
        R(j, j) = sum;
        // scale Q
        for (int i = 0; i < Q.rowsize(); i++) Q(i, j) = Q(i, j) / sum;
    }
    return R;
}

//! computes qr decomposition using modified Gram-Schmidt - works only with small (m>n) and square matrices
template <class T>
DenseMatrix<T> qr_gram_schmidt(DenseMatrix<T>& Q) {
    // create matrix R
    DenseMatrix<T> R(Q.colsize(), Q.colsize());

    // start orthogonalizing
    for (int k = 0; k < Q.colsize(); k++) {
        // compute norm of column k
        T sum_denom(0.0);
        for (int i = 0; i < Q.rowsize(); i++) {
            sum_denom += Q(i, k) * Q(i, k);
        }

        // fill the main diagonal of R with elements
        sum_denom = sqrt(sum_denom);
        R(k, k) = sum_denom;

        // scale column k to the main diagonal
        for (int i = 0; i < Q.rowsize(); i++) {
            Q(i, k) /= R(k, k);
        }

        // modify all later columns with column k
        for (int j = k + 1; j < Q.colsize(); j++) {
            // compute norm of column j
            T sum_nom(0.0);
            for (int i = 0; i < Q.rowsize(); i++) {
                sum_nom += Q(i, k) * Q(i, j);
            }
            // insert missing elements to R
            R(k, j) = sum_nom;

            // orthogonalize column j
            for (int i = 0; i < Q.rowsize(); i++) {
                Q(i, j) -= Q(i, k) * R(k, j);
            }
        }
    }
    return R;
}

//! computes qr decomposition using modified Gram-Schmidt and pivoting - works with all types of matrices
template <class T>
DenseMatrix<T> qr_gram_schmidt_pivoting(DenseMatrix<T>& Q, Vector<int>& p, int& rank, double threshold=0.00000000001) {
    // check if permutation vector has the right size
    if (p.size() != Q.colsize()) {
        HDNUM_ERROR("Permutation Vector incompatible with Matrix!");
    }

    // initialize permutation vector
    for (int i = 0; i < p.size(); i++) {
        p[i] = i;
    }

    // initialize rank
    rank = 0;

    // save matrix A, before it's replaced with Q
    DenseMatrix<T> A(Q);

    // create Matrix R
    hdnum::DenseMatrix<T> R(A.colsize(), A.colsize());

    // start orthogonalizing
    for (int k = 0; k < Q.colsize(); k++) {
        // find column with highest norm
        // compute norm of column k
        T norm_k(0.0);
        for (int r = 0; r < Q.rowsize(); r++) {
                norm_k += Q(r, k) * Q(r, k);
            }
        norm_k = sqrt(norm_k);

        // compare norm of column k to the following column norms
        for (int c = k+1; c < Q.colsize(); c++) {
            T norm(0.0);
            for (int r = 0; r < Q.rowsize(); r++) {
                norm += Q(r, c) * Q(r, c);
            }
            norm = sqrt(norm);
            // store permutation
            if (norm > norm_k) {
                p[k] = c;
            }
        }

        // swap columns if necessary
        if (p[k] > k) {
            for (int r = 0; r < Q.rowsize(); r++) {
                T temp_Q = Q(r, k);
                Q(r, k) = Q(r, p[k]);
                Q(r, p[k]) = temp_Q;
            }
            p[p[k]] = k;

            // compute norm of the new column k
            norm_k = 0;
            for (int i = 0; i < Q.rowsize(); i++) {
                norm_k += Q(i, k) * Q(i, k);
            }
            norm_k = sqrt(norm_k);
        }

        // if norm of column k > threshold -> column k is linear independent
        if (norm_k > threshold) {
            rank++;
        } else {
            break;
        }

        // modify all later columns with column k
        for (int j = k + 1; j < Q.colsize(); j++) {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i = 0; i < Q.rowsize(); i++) {
                sum_nom += Q(i, j) * Q(i, k);
                sum_denom += Q(i, k) * Q(i, k);
            }

            T alpha = sum_nom / sum_denom;
            for (int i = 0; i < Q.rowsize(); i++) Q(i, j) -= alpha * Q(i, k);
        }
    }

    // add values to R, except main diagonal
    for (int i = 1; i < R.colsize(); i++) {
        for (int j = 0; j < i; j++) {
            T sum_nom(0.0);
            T sum_l2nom(0.0);
            for (int k = 0; k < Q.rowsize(); k++) {
                sum_nom += A(k, p[i]) * Q(k, j);
                sum_l2nom += Q(k, j) * Q(k, j);
            }
            sum_l2nom = sqrt(sum_l2nom);
            // add element
            R(j, i) = sum_nom / sum_l2nom;
        }
    }

    // add missing values and scale
    for (int j = 0; j < Q.colsize(); j++) {
        // compute norm of column j
        T sum(0.0);
        for (int i = 0; i < Q.rowsize(); i++) sum += Q(i, j) * Q(i, j);
        sum = sqrt(sum);
        // add main diagonal to R
        R(j, j) = sum;
        // scale Q
        for (int i = 0; i < Q.rowsize(); i++) Q(i, j) = Q(i, j) / sum;
    }

    return R;
}

//! applies a permutation vector to a matrix
template<typename T>
void apply_permutation(DenseMatrix<T>& A, Vector<int>& p) {
    // check if permutation vector has the right size
    if (p.size() != A.colsize()) {
        HDNUM_ERROR("Permutation Vector incompatible with Matrix!");
    }

    // permutate the columns
    for (int k = 0; k < p.size(); k++) {
        if (p[k] != k) {
            // swap column
            for (int j=0; j < A.rowsize(); j++) {
                T temp_A = A(j, k);
                A(j, k) = A(j, p[k]);
                A(j, p[k]) = temp_A;
            }
            // swap inside permutation vector
            int temp_p = p[k];
            p[k] = p[temp_p];
            p[temp_p] = temp_p;
        }
    }
}

}  // namespace hdnum

#endif
