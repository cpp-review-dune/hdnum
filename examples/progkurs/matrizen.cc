// matrizen.cc
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

int main ()
{
  // Konstruktion
  hdnum::DenseMatrix<float> A(10,10);     // 10x10 Matrix uninitialisiert
  hdnum::DenseMatrix<float> B(10,10,1.0); // 10x10 Matrix uninitialisiert
  hdnum::DenseMatrix<float> C(10,10,2.0); // 10x10 Matrix initialisiert
  
  A[1][3] = 3.14; // Zuweisung auf A13
  A = B + C;      // Addition 
  C.mm(A,B);      // C = A*B (matmul)

  hdnum::Vector<float> x(10,1.0);
  hdnum::Vector<float> y(10,1.0);
  A.mv(y,x);    //y = A*x (matrix-vector)

  return 0;
}
