#include <iostream>
#include <vector>
#include "hdnum.hh"
#include <string>


using namespace hdnum;
 typedef double Number;               // define a number type

int main (){

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
  int s = 3;

  // compute invers of A (Ax_i=e_i  --> [x_1...x_s]=A⁻1)
  DenseMatrix<Number> Ainv (s,s,Number(0));
	for (int i=0; i < s; i++)                       
  {
     Vector<Number> e (s, Number(0));
     e[i]=Number(1);
     Vector<Number> w (s, Number(0));
     Vector<Number> x (s, Number(0));
     Vector<Number> y (s, Number(0));
     Vector<Number> z (s, Number(0));
     Vector<std::size_t> p(s);
     Vector<std::size_t> q(s);
DenseMatrix<Number> Temp (s,s,0.0);
Temp = A;
     row_equilibrate(A,w);                         // equilibrate rows
     lr_fullpivot(A,p,q);                          // LR decomposition of A
     apply_equilibrate(w,e);                       // equilibration of right hand side
     permute_forward(p,e);                         // permutation of right hand side
     solveL(A,e,e);                                // forward substitution
     solveR(A,z,e);                                // backward substitution
     permute_backward(q,z);                        // backward permutation
     for (int j = 0; j < s; j++)
     {
	     Ainv[i][j] = z[j];
     }
A = Temp;
  }
  std::cout << Ainv << std::endl; 

  return 0;
}

