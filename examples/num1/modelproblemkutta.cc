#include <iostream>
#include <vector>
#include "hdnum.hh"

#include "modelproblem.hh"
#include "../../src/testode.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef ModelProblem<Number> Model;  // Model type
  Model model(-1.0);                   // instantiate model


	/*hdnum::DenseMatrix<double> A(2,2,0.0);
	A[1][0] = 0.5;
 	A[1][1] = 0.5;

	hdnum::Vector<double> B(2,0.0);
	B[0] = 0.5;
	B[1] = 0.5;


	hdnum::Vector<double> C(3, 0.0);
	C[1] = 1; */
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

	
  typedef hdnum::RungeKutta_n<Model> Solver; // Solver type
  Solver solver(model, A, B, C);                // instantiate solver
  solver.set_dt(0.02);                  // set initial time step


  hdnum::Vector<Number> times;           // store time values here
  hdnum::Vector<hdnum::Vector<Number> > states; // store states here
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state
  while (solver.get_time()<5.0-1e-6) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
    }

  gnuplot("modelkutta.dat",times,states); // output model result

  return 0;
}
