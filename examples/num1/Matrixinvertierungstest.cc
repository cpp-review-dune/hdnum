#include <iostream>
#include <vector>
#include "hdnum.hh"
#include <string>


using namespace hdnum;

int main (){

  hdnum::DenseMatrix<double> A(2,2,0.0);
	A[1][0] = 1.0;
	A[0][1] = 1.0;

  int s = 2;

  // compute invers of A (Ax_i=e_i  --> [x_1...x_s]=A⁻1)
  DenseMatrix<double> Ainv (s,s,double(0));
	for (int i=0; i < s; i++)                       
  {
     Vector<double> e (s, double(0));
     e[i]=double(1);
     Vector<double> x (s, double(0));
     Vector<double> y (s, double(0));
     Vector<double> z (s, double(0));
     Vector<std::size_t> p(s);
     Vector<std::size_t> q(s);
std::cout << "Test1" << std::endl;
     row_equilibrate(A,y);                         // equilibrate rows
std::cout << "Test2" << std::endl;
     lr_fullpivot(A,p,q);                          // LR decomposition of A
std::cout << "Test3" << std::endl;
      apply_equilibrate(y,e);                       // equilibration of right hand side
std::cout << "Test4" << std::endl;
     permute_forward(p,e);                         // permutation of right hand side
std::cout << "Test5" << std::endl;
     solveL(A,e,e);                                // forward substitution
std::cout << "Test6" << std::endl;
     solveR(A,z,e);                                // backward substitution
std::cout << "Test7" << std::endl;
     permute_backward(q,z);                        // backward permutation
std::cout << "Test8" << std::endl;
     for (int j = 0; j < s; j++)
     {
	     Ainv[j][i] = z[j];
     }
 	}
  std::cout << Ainv << std::endl; 

  return 0;
}

