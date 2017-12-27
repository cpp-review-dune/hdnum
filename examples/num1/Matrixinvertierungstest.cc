#include <iostream>
#include <vector>
#include "hdnum.hh"
#include <string>


using namespace hdnum;

int main (){

  hdnum::DenseMatrix<double> A(3,3,0.0);
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
  DenseMatrix<double> Ainv (s,s,double(0));
	for (int i=0; i < s; i++)                       
  {
     Vector<double> e (s, double(0));
     e[i]=double(1);
     Vector<double> w (s, double(0));
     Vector<double> x (s, double(0));
     Vector<double> y (s, double(0));
     Vector<double> z (s, double(0));
     Vector<std::size_t> p(s);
     Vector<std::size_t> q(s);
DenseMatrix<double> Temp (s,s,0.0);
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

/*
  
  Vector<number> x(n);
  Vector<number> b(n);
  Vector<number> s(n);
  Vector<std::size_t> p(n);
  Vector<std::size_t> q(n);
  DenseMatrix<number> A(n,n);
  fill(x,number(1.0),number(1.0));
  vandermonde(A,x);
  A.mv(b,x);
  row_equilibrate(A,s);
  x = number(0.0);
  lr_fullpivot(A,p,q);
  apply_equilibrate(s,b);
  permute_forward(p,b);
  solveL(A,b,b);
  solveR(A,x,b);
  permute_backward(q,x);

*/

	
