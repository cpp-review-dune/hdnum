#include <iostream>
#include <vector>
#include "hdnum.hh"
#include <cmath>
#include "modelproblem.hh"


int main ()
{
  typedef double Number;               // define a number type
  typedef ModelProblem<Number> Model;  // Model type
  Model model(-1.0);                   // instantiate model
/*
  // implicit euler with s=1
  hdnum::DenseMatrix<double> A(1,1,0.0);
  A[0][0] = 1.0;
  
  hdnum::Vector<double> B(1,0.0);
  B[0] = 1.0;
  
  hdnum::Vector<double> C(1, 0.0);
  C[0] = 1.0; 
*/

  // Gauß with s=3
  hdnum::DenseMatrix<Number> A(3,3,0.0);
  A[0][0] = 5.0/36.0;
  A[0][1] = (10.0-3.0*sqrt(15.0))/45.0;
  A[0][2] = (25.0-6.0*sqrt(15.0))/180.0;
  A[1][0] = (10.0+3.0*sqrt(15.0))/72.0;
  A[1][1] = 2.0/9.0;
  A[1][2] = (10.0-3.0*sqrt(15.0))/72.0;
  A[2][0] = (25.0+6.0*sqrt(15.0))/180.0;
  A[2][1] = (10.0+3.0*sqrt(15.0))/45.0;
  A[2][2] = 5.0/36.0;
  
  hdnum::Vector<Number> B(3,0.0);
  B[0] = 5.0/18.0;
  B[1] = 4.0/9.0;
  B[2] = 5.0/18.0;

  hdnum::Vector<Number> C(3, 0.0);
  C[0] = (5.0-sqrt(15.0))/10.0;
  C[1] = 0.5;
  C[2] = (5.0+sqrt(15.0))/10.0;

  ordertest(model, A, B, C, 5.0, 0.2, 10);

  return 0;
}
