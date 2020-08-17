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
        for (int j=k+1; j<Q.rowsize(); j++)
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

}
#endif
