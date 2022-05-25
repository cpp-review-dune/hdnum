#ifndef NEWTONMETHODS_CC
#define NEWTONMETHODS_CC

#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include <cmath> 

#include "hdnum.hh"    // hdnum header
#include "nonlinearFunctions.hh" 
#include "newtonVisualization.hh" //contains test functions to visualize newton methods

using namespace hdnum;

int main ()
{
  // Examples of solving nonlinear models with the newton-raphson method:

  SquareRootProblem<double> sqp(16); // Squared root problem: min_x x^2 - a. (Here, a=16)
  std::cout<< "Size of the model: " << sqp.size()<<std::endl;

  // Declare a Newton-Raphson solver and use default valued for the maximum nuber of iterations, 
  // the number of line search steps, the absolute limit for defect and the reduction factor.
  Newton newtonRapshon;

  newtonRapshon.set_verbosity(1); // output the summary

  Vector<double> x(1);  // solution of the nonlinear problem will be saved here.
  x[0] = 25; // change the initial value.
  
  newtonRapshon.solve(sqp,  x);

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  Vector<double> solution(2); // Solution of the nonlinear problem will be saved here.

  // General nonlinear problems(For the loss function, you need to define a lambda function)
  auto problemRosenbrock = getNonlinearProblem(functionRosenbrock<double>,solution);
  auto problemBohachesky = getNonlinearProblem(functionBohachesky<double>, solution);
  auto problemBooth = getNonlinearProblem(functionBooth<double>, solution);
  auto problemBranin = getNonlinearProblem(functionBranin<double>, solution);
  auto problemMatyas = getNonlinearProblem(functionMatyas<double>, solution);

  std::cout<< "Number of components of model: " << problemRosenbrock.size()<<std::endl;

  // change settings of the solver
  newtonRapshon.set_linesearchsteps(50);
  newtonRapshon.set_maxit(10000);
  newtonRapshon.set_reduction(1e-15);

  solution[0] = -5;
  solution[1] = 2;

  newtonRapshon.solve(problemRosenbrock, solution);

  std::cout<< "Solution Newton: " << solution[0] << " " << solution[1] << std::endl;

  std::cout<<"-----------------------------------------------------------"<<std::endl;
  
  NewtonDogLegCauchy ndlc;
  solution[0] = -5;
  solution[1] = -2;
  ndlc.set_verbosity(1);
  ndlc.set_maxit(1000);
  ndlc.set_reduction(1e-15);
  ndlc.setInitialTrustRadius(100.0);
  ndlc.setMaxTrustRadius(100.0);
  ndlc.solve(problemRosenbrock,solution);
  std::cout<<"Solution: Newton Dog Cauchy: " << solution[0] << " " << solution[1] << std::endl;

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  // Viszualize newton and newton dogleg cauchy method
  ndlc.set_verbosity(0);
  newtonRapshon.set_verbosity(0);
  
  Vector<double> sol(2); // solutions of the solvers will be saved here
  sol[0] = -5.0;
  sol[1] = -2.0;
  
  Vector<double> domain = {-5.0,15.0};
  Vector<double> range = {0, 10};

  //testNewtonAgainstDogLeg(problemBranin, sol, domain, 10,10);
  //testDoglegConvergenceFixedRadius(problemBranin, 1.0,1.0, domain);
  //testDoglegConvergenceFixedInitialSolution(problemRosenbrock,sol, range);
  //testNewtonConvergence(problemBranin, domain);

  // Test Projected Newton method against different nonlinear minimization problems with constraints
  Vector<double> s(2);
  ProjectedNewton proj;

  DenseMatrix<double> constraints1(2,2);
  constraints1[0][0] = -1;
  constraints1[0][1] = 2;
  constraints1[1][0] = 1;
  constraints1[1][1] = 2;

  Vector<double> upperbound1 = {2,6};
  Vector<double> lowerbound1 = {-10000, -10000};
  s[0] = 0;
  s[1] = 1;
  auto prob1= getNonlinearMinimizationProblem_Constrained(functionConstrained1<double>, gradientConstrained1<double>, constraints1, lowerbound1, upperbound1, s);

  //proj.solve(prob1, s);

  DenseMatrix<double> constraints2(2,2);
  constraints2[0][0] = -2;
  constraints2[0][1] = -1;
  constraints2[1][0] = 1;
  constraints2[1][1] = 0;

  Vector<double> upperbound2 = {10000,10000};
  Vector<double> lowerbound2 = {-2, 0};
  s[0] = 0;
  s[1] = 1;
  auto prob2= getNonlinearMinimizationProblem_Constrained(functionConstrained2<double>, gradientConstraiend2<double>, constraints2, lowerbound2, upperbound2, s);

  //proj.solve(prob2, s);

  hdnum::DenseMatrix<double> constraints3(4,2);
  constraints3[0][0] = -1;
  constraints3[0][1] = 0;

  constraints3[1][0] = 0;
  constraints3[1][1] = -1;

  constraints3[2][0] = 1;
  constraints3[2][1] = 1;

  constraints3[3][0] = -1;
  constraints3[3][1] = 2;


  hdnum::Vector<double> upperbound3 = {0,0,8,10};
  hdnum::Vector<double> lowerbound3 = {-10000, -10000, -10000, -10000};
  s[0] = 0;
  s[1] = 0;
  auto prob3= hdnum::getNonlinearMinimizationProblem_Constrained(functionConstrained3<double>, gradientConstraiend3<double>, constraints3, lowerbound3, upperbound3, s);

  //proj.solve(prob3, s);

  //testProjectedNewton(prob3, constraints3, lowerbound3, upperbound3);
}

#endif


