#include <iostream>
#include <vector>
#include "hdnum.hh"
#include "modelproblem.hh"


int main ()
{
  typedef double Number;               // define a number type

  typedef ModelProblem<Number> Model;  // Model type
  Model model(-1.0);                   // instantiate model

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

  typedef hdnum::RungeKutta_n<Model> Solver;    // Solver type
  Solver solver(model, A, B, C);                // instantiate solver
  solver.set_dt(Number(0.02));                  // set initial time step


  hdnum::Vector<Number> times;                  // store time values here
  hdnum::Vector<hdnum::Vector<Number> > states; // store states here
  times.push_back(solver.get_time());           // initial time
  states.push_back(solver.get_state());         // initial state
  while (solver.get_time()<1.04-1e-6)            // the time loop
  {
    solver.step();                              // advance model by one time step
    times.push_back(solver.get_time());         // save time
    states.push_back(solver.get_state());       // and state
  }

  gnuplot("modelkutta.dat",times,states); // output model result

  return 0;
}
