#ifndef NEWTONMETHODS_CC
#define NEWTONMETHODS_CC

#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include <cmath> 

#include "hdnum.hh"    // hdnum header* rightPart 
#include "nonlinearFunctions.hh" 
#include "newtonVisualization.hh" //contains test functions to visualize newton methods
#include <any>
#include <algorithm>

using namespace hdnum;


int main ()
{
  // Examples of solving nonlinear models with the newton-raphson method:

  std::cout<<"Example of solving squared root problem with newton method: \n "<<std::endl;
  SquareRootProblem<double> sqp(16); // Squared root problem: min_x x^2 - a. (Here, a=16)

  // Declare a Newton-Raphson solver and use default valued for the maximum nuber of iterations, 
  // the number of line search steps, the absolute limit for defect and the reduction factor.
  Newton newtonRapshon;

  newtonRapshon.set_verbosity(1); // output the summary

  Vector<double> x{25};  // solution of the nonlinear problem will be saved here.
  
  newtonRapshon.solve(sqp,  x);

  std::cout<<"\nExample of solving general nonlinear problems with newton and newton dog leg cauchy method: \n" << std::endl;
  Vector<double> solution(2); // Solution of the nonlinear problem will be saved here.

  // General nonlinear problems(For the loss function, you need to define a lambda function)
  auto problemRosenbrock = getNonlinearProblem(&functionRosenbrock<double>,solution);
  auto problemBohachesky = getNonlinearProblem(&functionBohachesky<double>, solution);
  auto problemBooth = getNonlinearProblem(&functionBooth<double>, solution);
  auto problemBranin = getNonlinearProblem(&functionBranin<double>, solution);
  auto problemMatyas = getNonlinearProblem(&functionMatyas<double>, solution);

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
  ndlc.set_maxit(20);
  ndlc.set_abslimit(1e-14);
  ndlc.set_reduction(1e-15);
  ndlc.setMaxTrustRadius(100.0);
  ndlc.setInitialTrustRadius(100.0);
  ndlc.solve(problemRosenbrock,solution);
  std::cout<<"Solution: Newton Dog Cauchy: " << solution[0] << " " << solution[1] << std::endl;

  std::cout<<"\nExample of solving general complex nonlinear problems with the newton and newton dog leg cauchy method: \n" << std::endl;

  CVector<double> complexSolution(1);
  complexSolution[0].real(8.0);
  complexSolution[0].imag(4.0);

  auto complexproblem = getNonlinearProblem(&complexFunction<double>, complexSolution);
  newtonRapshon.solve(complexproblem, complexSolution);
  std::cout<<"Solution: " << complexSolution[0] <<std::endl;

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  complexSolution[0].real(8);
  complexSolution[0].imag(4);
  ndlc.solve(complexproblem, complexSolution);
  std::cout<<"complexSolution: " << complexSolution[0] <<std::endl;
  
  // Viszualize newton and newton dogleg cauchy method
  newtonRapshon.set_verbosity(0);
  ndlc.set_verbosity(0);
  
  Vector<double> sol(2); // solutions of the solvers will be saved here
  sol[0] = -5.0;
  sol[1] = -2.0;
  
  Vector<double> domain = {-5.0,15.0};
  Vector<double> range = {0, 10};

  // testNewtonAgainstDogLeg(problemRosenbrock, sol, domain, 10,10);
  // testDoglegConvergenceFixedRadius(problemBranin, 1.0,1.0, domain);
  // testDoglegConvergenceFixedInitialSolution(problemRosenbrock,sol, range);
  // testNewtonConvergence(problemBranin, domain);

  Vector<std::complex<double>> complexsol(1);
  complexsol[0].real(8.0);
  complexsol[0].imag(4.0);
  //testNewtonAgainstDogLegLossAndReduction(complexproblem, complexsol, 1.0,1.0);

  // Examples of solving constrained nonlinear minimization problem with the projected Newton method
  Vector<double> projectedNewtonSolution(2);
  projectedNewtonSolution[0] = 0;
  projectedNewtonSolution[1] = 1;

  auto constraintProblem1= getNonlinearMinimizationProblem_Constrained(&ConstrainedProblem1<double>::objective, &ConstrainedProblem1<double>::gradient, &ConstrainedProblem1<double>::hessian, ConstrainedProblem1<double>::constraints(), ConstrainedProblem1<double>::lowerbounds(), ConstrainedProblem1<double>::upperbounds(), projectedNewtonSolution);
  auto constraintProblem2= getNonlinearMinimizationProblem_Constrained(&ConstrainedProblem2<double>::objective, &ConstrainedProblem2<double>::gradient, &ConstrainedProblem2<double>::hessian, ConstrainedProblem2<double>::constraints(), ConstrainedProblem2<double>::lowerbounds(), ConstrainedProblem2<double>::upperbounds(), projectedNewtonSolution);
  auto constraintProblem3= getNonlinearMinimizationProblem_Constrained(&ConstrainedProblem3<double>::objective, &ConstrainedProblem3<double>::gradient, &ConstrainedProblem3<double>::hessian, ConstrainedProblem3<double>::constraints(), ConstrainedProblem3<double>::lowerbounds(), ConstrainedProblem3<double>::upperbounds(), projectedNewtonSolution);

  Vector<double> projectedNewtonDomain = {-5,10};
  Vector<double> projectedNewtonInitialSolution = {-2,6};
 
  // testProjectedNewton(constraintProblem1, projectedNewtonInitialSolution, constraints1, lowerbound1, upperbound1, projectedNewtonDomain);
  // testProjectedNewton(constraintProblem2, projectedNewtonInitialSolution, constraints2, lowerbound2, upperbound2, projectedNewtonDomain);
  // testProjectedNewton(constraintProblem3, projectedNewtonInitialSolution, constraints3, lowerbound3, upperbound3, projectedNewtonDomain);

  std::cout <<"\nMinimize distance to circle on the domain in a shape of an infinite pyramid" << std::endl;

  Vector<double> minDistToCirlceSolution = {12, 12, 20};
  auto minDistToCirlceProblem = getNonlinearMinimizationProblem_Constrained(&MinDistToCircle<double>::objective, &MinDistToCircle<double>::gradient , &MinDistToCircle<double>::hessian , MinDistToCircle<double>::constraints(), MinDistToCircle<double>::lowerbounds(), MinDistToCircle<double>::upperbounds(), minDistToCirlceSolution);

  ProjectedNewton projnewt;
  projnewt.set_maxit(20);
  projnewt.solve(minDistToCirlceProblem, minDistToCirlceSolution);
  std::cout<< "Solution: "<< minDistToCirlceSolution<< std::endl;
}
#endif


