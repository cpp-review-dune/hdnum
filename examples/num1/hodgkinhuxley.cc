#include <iostream>
#include <vector>
#include "hdnum.hh"

using namespace hdnum;

#include "hodgkinhuxley.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef HodgkinHuxley<Number> Model; // Model type
  Model model;                         // instantiate model

  Number TOL=1e-4; // desired tolerance
  
  // adaptive embedded RK method
  typedef RKF45<Model> Solver;         // Solver type
  Solver solver(model);                // instantiate solver
  solver.set_dt(1.0/16.0);             // set initial time step
  solver.set_TOL(TOL);

  // Richardson Extrapolation with RK4
  // typedef RungeKutta4<Model> SubSolver;         // Solver type
  // //typedef Heun2<Model> SubSolver;         // Solver type
  // SubSolver subsolver(model);                // instantiate solver
  // typedef RE<Model,SubSolver> Solver;         // Solver type
  // Solver solver(model,subsolver);                // instantiate solver
  // solver.set_dt(1.0/16.0);             // set initial time step
  // solver.set_TOL(TOL); 
  // solver.set_dtmin(1e-8); 
  // solver.set_rho(0.8); 

  std::vector<Number> times;           // store time values here
  std::vector<Vector<Number> > states; // store states here
  std::vector<Number> dts;             // store delta t
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state
  dts.push_back(solver.get_dt());      // initial dt

  Number T = 130.0;
  while (solver.get_time()<T-1e-6) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
      dts.push_back(solver.get_dt());      // used dt
      //std::cout << solver.get_time() << " " << solver.get_dt() << std::endl;
    }

  gnuplot("hodgkinhuxley_rk45.dat",times,states,dts); // output model result

  std::cout << times.size() << " timesteps" << std::endl;
  std::cout << model.get_count() << " f evaluations" << std::endl;
  
  return 0;
}
