#include <iostream>
#include <vector>
#include "hdnum.hh"
#include "modelproblem.hh"


int main ()
{
  typedef double Number;               // define a number type

  typedef ModelProblem<Number> Model;  // Model type
  Model model(-1.0);                   // instantiate model

  // implicit solvers

/*
	// implicit euler
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

/*
	//implicit Alexander Dirk Verfahren
  double alpha = 1+sqrt(2)/2;
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

  // Gauß with s=3
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

	hdnum::Vector<double> B(3,0.0);
	B[0] = 5.0/18.0;
	B[1] = 4.0/9.0;
	B[2] = 5.0/18.0;

	hdnum::Vector<double> C(3, 0.0);
	C[0] = (5.0-sqrt(15.0))/10.0;
	C[1] = 0.5;
	C[2] = (5.0+sqrt(15.0))/10.0;


	// explicit solver
	
/*	
  // classical runge kutta method
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
*/	
	
  typedef hdnum::RungeKutta_n<Model> Solver;    // Solver type
  Solver solver(model, A, B, C);                // instantiate solver
  solver.set_dt(Number(0.02));                  // set initial time step


  hdnum::Vector<Number> times;                  // store time values here
  hdnum::Vector<hdnum::Vector<Number> > states; // store states here
  times.push_back(solver.get_time());           // initial time
  states.push_back(solver.get_state());         // initial state
  while (solver.get_time()<5.0-1e-6)            // the time loop
  {
    solver.step();                              // advance model by one time step
    times.push_back(solver.get_time());         // save time
    states.push_back(solver.get_state());       // and state
  }

  gnuplot("modelkutta.dat",times,states); // output model result

  return 0;
}
