#ifndef HDNUM_QRHOUSHOLDER_HH
#define HDNUM_QRHOUSHOLDER_HH
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include "vector.hh"
#include "densematrix.hh"
#include <cmath>
/** @file
 *  @brief This file implements QR decomposition using housholder transformation
 */
namespace hdnum {
    template<class T>
    Vector<T> Calculate_Ai(DenseMatrix<T>& A){  
              Vector<T> a(A.rowsize(),0);
        for (int i =0;i<A.rowsize();i++){
            a[i]=A(i,0);
        }
        T alpha= sqrt(a.two_norm_2());
        a[0]=a[0]-alpha;

        return a;
    }
          /*! \brief Function that calculate the Vi vector of the firss column of matrix A
           */
    template<class T>
    Vector<T> Vi(DenseMatrix<T>& A){
        Vector<T> a(A.rowsize(),0);
        for (int i =0;i<A.rowsize();i++){
            a[i]=A(i,0);
        }
        T alpha= sqrt(a.two_norm_2());
        a[0]=a[0]-alpha;
        return a;
    }
              /*! \brief Function that gives the sign of a number
           */
    template <typename REAL> int sgn(REAL val) {
    return (REAL(0) < val) - (val < REAL(0));
}
      /*! \brief Funktion that calculate the QR decoposition in place 
      the elements of A will be replaced with the elements of R and the under diagonals elements of
      vectors V1-V2...Vi and the diagonal elements of vectors Vi will be saved in vectro B.
      */
     template<class REAL>
  void qrhousholder (DenseMatrix<REAL>& A, std::vector<REAL>& v)
  {
    auto m = A.rowsize();
    auto n = A.colsize();
    DenseMatrix<REAL> u(m,n,0);
    if(n>m){
      std::cout<<"In the QR decompesition Rang(A) = rowssize but here colsize > rowsize try another matrix";
      return;
    }
    for (std::size_t k=0; k<n ; k++){
      REAL alpha = 0;
      REAL a = abs(A(k,k));
      // calculate the H_{i} matrix
      for(std::size_t i= k;i<m;i++){
        u(i,k)=A(i,k)/a;
        alpha = alpha + pow(u(i,k),2);
      }
      alpha = sqrt(alpha);
      REAL B = 1/(alpha * (alpha+abs(u(k,k))));
      u(k,k)=u(k,k)+(sgn(A(k,k)) * alpha);
      //multiplicate H_{k} with A
      A(k,k)=-sgn(A(k,k)) * a * alpha;
      for (std::size_t i= k+1;i<=m;i++){
        A(i,k)=0;
      }
      for (std::size_t j=k+1;j<=n;j++){
        REAL sum =0;
        for(std::size_t i=k;i<=m;i++){
          sum=sum + (u(i,k) * A(i,k)); 
        }
        REAL s= B * sum;
        for(std::size_t i = k ; i <= m;i++){
          A(i,j)=A(i,j)-(s * u(i,k));
        }
      }
    }
  }
     template<class REAL>
  void qrhousholder1 (DenseMatrix<REAL>& A, std::vector<REAL>& v)
  {
    auto m = A.rowsize();
    auto n = A.colsize();
    //std::vector<REAL>temp(n,0);
    //for (int i=0;i<n;i++){
     // temp[i]=A(i,i);
    //}
     for (int j=0;j<n;j++){
       REAL s=0;
       for(int i=j;i<m;i++){
         s=s+pow(A(i,j),2);
       }
       s=sqrt(s);
       v[j] = (-1.0) * sgn(A(j,j)) * s;
       REAL fak = sqrt( s * ( s + std::abs(A(j,j))));
       A(j,j)=A(j,j)-v[j];
       for (int k=j ;k< m;k++){
         A(k,j)=A(k,j)/fak;
       }
       for (int i =j+1 ;i<n;i++){
         s=0;
         for (int k=j;k<m;k++){
           s=s+ (A(k,j) * A(k,i));
         }
         for (int k=j;k<m;k++){
           A(k,i)=A(k,i)- ( A(k,j) * s );
         }
       }
       //normalize the vi vectors again
       for (int i=m;i>=0;i--){
         A(i,j)=A(i,j)*fak;
         if(i==j){
           break;
         }
       }
     }
     //for (int i=0;i<n;i++){
     //  for(int j=m;j>0;j--){
     //    A(j,i)=A(j,i) * temp[i];
     //    if(i==j){
       //    break;
        // }
       //}
     //}
  }
  template<class REAL>
  DenseMatrix<REAL>qrhousholderexplizitQ (DenseMatrix<REAL>& A, std::vector<REAL>& v)
  {
    auto m = A.rowsize();
    auto n = A.colsize();
    //std::vector<REAL>temp(n,0);
    //for (int i=0;i<n;i++){
     // temp[i]=A(i,i);
    //}
     for (int j=0;j<n;j++){
       REAL s=0;
       for(int i=j;i<m;i++){
         s=s+pow(A(i,j),2);
       }
       s=sqrt(s);
       v[j] = (-1.0) * sgn(A(j,j)) * s;
       REAL fak = sqrt( s * ( s + std::abs(A(j,j))));
       A(j,j)=A(j,j)-v[j];
       for (int k=j ;k< m;k++){
         A(k,j)=A(k,j)/fak;
       }
       for (int i =j+1 ;i<n;i++){
         s=0;
         for (int k=j;k<m;k++){
           s=s+ (A(k,j) * A(k,i));
         }
         for (int k=j;k<m;k++){
           A(k,i)=A(k,i)- ( A(k,j) * s );
         }
       }
       //normalize the vi vectors again
       for (int i=m;i>=0;i--){
         A(i,j)=A(i,j)*fak;
         if(i==j){
           break;
         }
       }
     }
    //q berechnen
     /* REAL y=0;
     for (int j=n ;j >=0;j++){
       REAL s= 0;
       for (int k=j;k<m;k++){
        if(k==j){
          y=1;
        }
        else if(k!=j){
          y=0;
        }
        s=s+ (A(k,j) * y);
       }
     } */
  }
}
#endif 