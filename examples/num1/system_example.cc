#include <iostream>
#include <vector>
#include "hdnum.hh"

#include "system_example.hh"

int main ()
{
  // Define a number type
  typedef double Number;

  // Model problem
  typedef SystemExample<Number> Model;
  Model model;

  // ODE solver
  typedef hdnum::EE<Model> Solver;
  Solver solver(model);
  solver.set_dt(0.02);

  // Store times and states
  hdnum::Vector<Number> times;
  hdnum::Vector<hdnum::Vector<Number>> states;
  times.push_back(solver.get_time());
  states.push_back(solver.get_state());

  // Time loop
  while (solver.get_time()<5.0-1e-6)
    {
      solver.step();
      times.push_back(solver.get_time());
      states.push_back(solver.get_state());
    }

  // Output
  gnuplot("system_example.dat",times,states);

  return 0;
}
