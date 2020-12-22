// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_QR_HH
#define HDNUM_QR_HH

#include "vector.hh"
#include "densematrix.hh"
#include <cmath>

/** @file
 *  @brief This file implements QR decomposition
 */

namespace hdnum
{
  //! computes orthonormal basis of Im(A) using classical Gram-Schmidt
  template<class T>
  DenseMatrix<T> gram_schmidt (const DenseMatrix<T>& A)
  {
    DenseMatrix<T> Q(A);

    // for all columns except the first
    for (int k=1; k<Q.colsize(); k++)
      {
        // orthogonalize column k against all previous
        for (int j=0; j<k; j++)
          {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i=0; i<Q.rowsize(); i++)
              {
                sum_nom += A[i][k]*Q[i][j];
                sum_denom += Q[i][j]*Q[i][j];
              }
            // modify
            T alpha = sum_nom/sum_denom;
            for (int i=0; i<Q.rowsize(); i++)
              Q[i][k] -= alpha*Q[i][j];
          }
      }
    for (int j=0; j<Q.colsize(); j++)
      {
        // compute norm of column j
        T sum(0.0);
        for (int i=0; i<Q.rowsize(); i++) sum += Q[i][j]*Q[i][j];
        sum = sqrt(sum);
        //scale
        for (int i=0; i<Q.rowsize(); i++) Q[i][j] = Q[i][j]/sum;
      }
    return Q;
  }

  //! computes orthonormal basis of Im(A) using modified Gram-Schmidt
  template<class T>
  DenseMatrix<T> modified_gram_schmidt (const DenseMatrix<T>& A)
  {
    DenseMatrix<T> Q(A);

    for (int k=0; k<Q.colsize(); k++)
      {
        // modify all later columns with column k
        for (int j=k+1; j<Q.colsize(); j++)
          {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i=0; i<Q.rowsize(); i++)
              {
                sum_nom += Q[i][j]*Q[i][k];
                sum_denom += Q[i][k]*Q[i][k];
              }
            // modify
            T alpha = sum_nom/sum_denom;
            for (int i=0; i<Q.rowsize(); i++)
              Q[i][j] -= alpha*Q[i][k];
          }
      }
    for (int j=0; j<Q.colsize(); j++)
      {
        // compute norm of column j
        T sum(0.0);
        for (int i=0; i<Q.rowsize(); i++) sum += Q[i][j]*Q[i][j];
        sum = sqrt(sum);
        //scale
        for (int i=0; i<Q.rowsize(); i++) Q[i][j] = Q[i][j]/sum;
      }
    return Q;
  }

  //! computes qr decomposition using modified Gram-Schmidt
  template<class T>
  DenseMatrix<T> qr_decomposition_gram_schmidt (DenseMatrix<T>& Q)
  {
    // check preconditions
    // Is A a small matrix?
    if (Q.colsize() >= Q.rowsize())
      {
        HDNUM_ERROR("The Matrix is not a small matrix!");
      }

    // save matrix A, before it's replaced with Q
    DenseMatrix<T> A(Q);

    // create matrix R
    DenseMatrix<T> R(Q.colsize(), Q.colsize());

    for (int k=0; k<Q.colsize(); k++)
      {
        // modify all later columns with column k
        for (int j=k+1; j<Q.colsize(); j++)
          {
            // compute factor
            T sum_nom(0.0);
            T sum_denom(0.0);
            for (int i=0; i<Q.rowsize(); i++)
              {
                sum_nom += Q(i, j)*Q(i, k);
                sum_denom += Q(i, k)*Q(i, k);
              }

            // add first value (1, 1) to R
            if (k == 0) R(0, 0) = sqrt(sum_denom);

            T alpha = sum_nom/sum_denom;
            for (int i=0; i<Q.rowsize(); i++)
              Q(i, j) -= alpha*Q(i, k);
          }
      }

    // add values to R, except main diagonal
    for (int i=1; i<R.colsize(); i++)
      {
        for (int j=0; j<i; j++)
          {
            T sum_nom(0.0);
            T sum_l2nom(0.0);
            for (int k=0; k<Q.rowsize(); k++)
              {
                sum_nom += A(k, i)*Q(k, j);
                sum_l2nom += Q(k, j)*Q(k, j);
              }
            sum_l2nom = sqrt(sum_l2nom);
            // add element
            R(j, i) = sum_nom/sum_l2nom;
          }
      }

    for (int j=0; j<Q.colsize(); j++)
      {
        // compute norm of column j
        T sum(0.0);
        for (int i=0; i<Q.rowsize(); i++) sum += Q(i, j)*Q(i, j);
        sum = sqrt(sum);
        // add main diagonal to R
        if (j > 0) R(j, j) = sum;
        // scale Q
        for (int i=0; i<Q.rowsize(); i++) Q(i, j) = Q(i, j)/sum;
      }

    return R;
  }
}

#endif
