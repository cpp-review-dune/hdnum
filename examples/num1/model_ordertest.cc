#include <iostream>
#include <vector>
#include "hdnum.hh"

#include "modelproblem.hh"


int main ()
{
  typedef double Number;               // define a number type

  typedef ModelProblem<Number> Model;  // Model type
  Model model(-1.0);                   // instantiate model

/*
	//implicit euler
	hdnum::DenseMatrix<double> A(1,1,0.0);
 	A[0][0] = 1.0;

	hdnum::Vector<double> B(1,0.0);
	B[0] = 1.0;


	hdnum::Vector<double> C(1, 0.0);
	C[0] = 1.0; 
*/

/*
	//implicit Trapez Regel
	hdnum::DenseMatrix<double> A(2,2,0.0);
	A[1][0] = 0.5;
 	A[1][1] = 0.5;

	hdnum::Vector<double> B(2,0.0);
	B[0] = 0.5;
	B[1] = 0.5;


	hdnum::Vector<double> C(2, 0.0);
	C[1] = 1; 
*/

	//implicit Alexander Dirk Verfahren
  /*	double alpha = 1+sqrt(2)/2;
	hdnum::DenseMatrix<double> A(2,2,0.0);
 	 A[0][0] = alpha;
	A[1][0] = 1-alpha;
 	A[1][1] = alpha;

	hdnum::Vector<double> B(2,0.0);
	B[0] = 1-alpha;
	B[1] = alpha;


	hdnum::Vector<double> C(2, 0.0);
  	C[0] = alpha;
	C[1] = 1; 
*/

	//explicit
	
	hdnum::DenseMatrix<double> A(4,4,0.0);
	A[1][0] = 0.5;
	A[2][1] = 0.5;
	A[3][2] = 1.0;

	hdnum::Vector<double> B(4,0.0);
	B[0] = 1.0/6.0;
	B[1] = 2.0/6.0;
	B[2] = 2.0/6.0;
	B[3] = 1.0/6.0;

	hdnum::Vector<double> C(4, 0.0);
	C[1] = 0.5;
	C[2] = 0.5;
	C[3] = 1.0;
	


  ordertest(model, A, B, C, 0.2, 0.2, 20);





  return 0;
}
