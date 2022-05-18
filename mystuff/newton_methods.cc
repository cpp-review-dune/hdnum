#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include "hdnum.hh"    // hdnum header
#include <vector>
#include <cmath>     

using namespace hdnum;

// arange types
using Number = double;
using Vec = Vector<Number>;
using Mat = DenseMatrix<Number>;

// Functions F(x) which describe the problem F(x) = 0
auto functionMatyas = [](const Vec& x) {
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = 2 * 0.26 * x[0] - 0.48 * x[1];
  result[1] = 2 * 0.26 * x[1] - 0.48 * x[0];

  return result;
};

auto functionRosenbrock = [](const Vec& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = -2.0 * (1-x[0]) - 400.0 * (x[1] - x[0] * x[0]) * x[0];
  result[1] = 200.0 * (x[1] -x[0] * x[0]);

  return result;
};

auto functionBohachesky = [](const Vec& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = 2.0 * x[0] + 0.3 * 3.0 * M_PI * sin(3.0 * M_PI * x[0]);
  result[1] = 4.0 * x[1] + 0.4 * 4.0 * M_PI * sin(4.0 * M_PI * x[1]);
  
  return result;
};

auto functionBooth = [](const Vec& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = 2.0 * (x[0] + 2.0 * x[1] - 7.0) + 4.0 * (2.0 * x[0] +x[1]- 5.0);
  result[1] = 4.0 * (x[0] + 2.0 * x[1] - 7.0) + 2.0 * (2.0 * x[0] +x[1]- 5.0);
  return result;
};
  
auto functionBranin = [](const Vec& x) 
{ 
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  Number a = 1.0;
  Number b = 5.1 / ( 4.0 * M_PI * M_PI );
  Number c = 5.0 / M_PI;

  Number r = 6.0;
  Number s = 10.0;
  Number t = 1.0 / (8 * M_PI);

  result[0] = 2.0 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r) * (-2.0 * b * x[0] + c) - s * (1 - t) * sin(x[0]);
  result[1] = 2.0 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r);

  return result;
};

//Functions F(x) which describe the problem min F(x) (min 0.5*||F(x)||^2 for higher dimensional features) 
auto functionConstrained1 = [](const Vec& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(1);
  result[0] = (x[0] - 1.0 ) * (x[0] - 1.0 ) + (x[1] - 2.5 ) *  (x[1] - 2.5 );

  return result;
};

auto gradientConstrained1 = [](const Vec& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = 2.0 * (x[0] - 1.0 );
  result[1] = 2.0 * (x[1] - 2.5 );

  return result;
};

auto functionConstrained2 = [](const Vec& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(1);
  result[0] = x[0]*x[0] + x[1]*x[1] - 4.0 * x[0] - 5.0 * x[1] + 2.0;

  return result;
};

auto gradientConstraiend2 = [](const Vec& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = 2.0 * x[0] - 4.0;
  result[1] = 2.0 * x[1] - 5.0;

  return result;
};

auto functionConstrained3=[](const Vec& x){
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(1);
  result[0] = 2.0 * x[0]*x[0] + x[1]*x[1] - 2.0 * x[0] * x[1] - 4.0 * x[0] - 6.0 * x[1];

  return result;
};

auto gradientConstraiend3 = [](const Vec& x)
{
  if(x.size() != 2){
    HDNUM_ERROR("Size of vector not equal 2");
  }
  Vec result(2);
  result[0] = 4.0 * x[0] - 2.0 * x[1] - 4.0;
  result[1] = 2.0 * x[1] - 2.0 * x[0] - 6.0;

  return result;
};


/**
 * @brief Test Newton method against Newton dogleg cauchy method(with a fixed maximum trust radius and initial trust radius).
 *        This test method will only work for the input dimension 2.
 *        Both methods output a summary on the terminal and afterwards some results will be vizualized:
 *        - trajectory of the optimal solution(3D and 2D)
 *        - process of norm and reduction
 *        - vizualizing the trust radii of the Newton dog leg cauchy method
 * 
 * @tparam NonlinearProblem 
 * @param nonlinearProblem 
 * @param initialSolution 
 * @param domain // domain for the plot
 * @param maxTrustRadius // default: 1.0
 * @param initialTrustRadius // default: 1.0
 */
template<typename NonlinearProblem>
void testNewtonAgainstDogLeg(const NonlinearProblem& nonlinearProblem, const Vector<typename NonlinearProblem::number_type>& initialSolution, const Vector<typename NonlinearProblem::number_type> domain, const double& maxTrustRadius= 1.0,const double& initialTrustRadius= 1.0){
    if(nonlinearProblem.size() != 2){
      HDNUM_ERROR("Error: This test method will only work for the input dimension 2");
    }
    typedef typename NonlinearProblem::number_type N;

    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<"Test Newton method against Newton Dogleg Cauchy"<< std::endl;
    std::cout<<"Initial solution: " << initialSolution[0] << " " << initialSolution[1]<< std::endl;
    std::cout<<"Domain region: " << domain[0] << " " << domain[1] << std::endl;
    std::cout<<"Initial trust radius: " << initialTrustRadius << " , Maxmimum trust radius: " << maxTrustRadius << std::endl;
    Vector<N> solution = initialSolution;
    
    std::fstream file_loss("loss.dat",std::ios::out);
    file_loss << domain[0]<<" " << domain[1] << "\n";

    N stepSize = 0.2;

    for(N i = domain[0]; i<= domain[1] ; i = i+stepSize){
      for(N j = domain[0]; j<= domain[1] ; j = j+stepSize){
        Vector<N> point = {i,j};
        Vector<N> result(2);
        nonlinearProblem.F(point, result);
        N value = 0.5 * (result * result); 
        file_loss<< value << "\n";
      }
    }
    file_loss.close();

    Newton newtonSolver;
    newtonSolver.set_verbosity(1);
    newtonSolver.set_linesearchsteps(50);
    newtonSolver.set_maxit(500);

    newtonSolver.solve(nonlinearProblem, solution, "newton_solver.dat");

    solution = initialSolution;

    NewtonDogLegCauchy doglegSolver;
    doglegSolver.set_verbosity(1);
    doglegSolver.set_maxit(500);
    doglegSolver.set_reduction(1e-15);
    doglegSolver.set_initial_trust_radius(initialTrustRadius);
    doglegSolver.set_max_radius(maxTrustRadius);

    doglegSolver.solve(nonlinearProblem, solution,"dog_leg_solver.dat");

    system("python3 newton_methods_test_against.py");

    std::remove("run.py");
    std::remove("loss.dat");
    std::remove("newton_solver.dat");
    std::remove("dog_leg_solver.dat");

    std::cout<<"-------------------------------"<<std::endl;
}

/**
 * @brief Test the convergence of the Newton dog leg cauchy method with different initial solutions and 
 *        a fixed initial and maximum trust radius. This test method will only work for the input dimension 2.
 * 
 * @tparam NonlinearProblem 
 * @tparam N number type   
 * @param initialTrustRadius 
 * @param maxTrustRadius 
 * @param nonlinearProblem 
 * @param domain 
 */
template<typename NonlinearProblem>
void testDoglegConvergenceFixedRadius(const NonlinearProblem& nonlinearProblem, const double& initialTrustRadius,const double& maxTrustRadius, Vector<typename NonlinearProblem::number_type> domain){
  if(nonlinearProblem.size() != 2){
    HDNUM_ERROR("Error: This test method will only work for the input dimension 2");
  }
  typedef typename NonlinearProblem::number_type N;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"Test Newton dog leg cauchy method for different initial solutions and a fixed initial and maximum trust radius"<< std::endl;
  std::cout<<"Initial trust radius: " << initialTrustRadius << " , maximum trust radius: " << maxTrustRadius << std::endl;
  std::cout<<"Domain: " << domain[0] << " " << domain[1] << std::endl;

  Vector<N> solution(nonlinearProblem.size());
  NewtonDogLegCauchy doglegSolver;
  doglegSolver.set_verbosity(0);
  doglegSolver.set_maxit(1000);
  doglegSolver.set_max_radius(maxTrustRadius);
  doglegSolver.set_initial_trust_radius(initialTrustRadius);

  std::fstream fileConvergence("convergence.dat",std::ios::out);
  fileConvergence <<"1" << "\n"; 
  fileConvergence << domain[0] << " " << domain[1]<<"\n";

  //change the step size to quick up process
  N step_size = 0.2;

  for(N i=domain[0]; i<=domain[1]; i=i+step_size){
    for(N j=domain[0]; j<=domain[1]; j=j+step_size){
      solution[0] = i;
      solution[1] = j;

      doglegSolver.solve(nonlinearProblem, solution);
      if(doglegSolver.has_converged())
        fileConvergence << doglegSolver.iterations()<< "\n";
      else
        fileConvergence << -1<< "\n";
    }
  }
  
  fileConvergence.close();

  system("python3 test_convergence.py");
  std::remove("convergence.dat");
  std::cout<<"-------------------------------"<<std::endl;
}

/**
 * @brief Test the convergence of the Newton dog leg cauchy method for different initial and maximum trust radii 
 *        and a fixed initial solution. This method will work for any input dimension.
 * 
 * @tparam NonlinearProblem 
 * @param initialSolution 
 * @param nonlinearProblem 
 * @param range // [minimum max_trust_radius, maximum max_trust_radius]
 */
template<typename NonlinearProblem>
void testDoglegConvergenceFixedInitialSolution(const NonlinearProblem& nonlinearProblem, const Vector<typename NonlinearProblem::number_type>& initialSolution, Vector<double> range){
  typedef typename NonlinearProblem::number_type N;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"Test Newton dog leg cauchy method for different initial and maximum trust radii and a fixed initial solution"<< std::endl;
  std::cout<<"Initial solution: ";
  for(auto sol: initialSolution){
    std::cout<< sol;
  }
  std::cout<<std::endl;
  std::cout<<"Range: " << range[0] << " " << range[1] << std::endl;

  Vector<N> solution(2);
  NewtonDogLegCauchy doglegSolver;

  std::fstream fileConvergence("convergence.dat",std::ios::out);
  fileConvergence <<"2" << "\n"; 
  fileConvergence << range[0] << " " << range[1]<<"\n";

  //change the step size to quick up process
  double stepSize = 0.2;

  for(double i=range[0]; i<=range[1]; i=i+stepSize){
    for(double j=range[0]; j<=range[1]; j=j+stepSize){
      if(j > i){
        fileConvergence << -1<< "\n";
        continue;
      }
      solution = initialSolution;

      doglegSolver.set_verbosity(0);
      doglegSolver.set_maxit(500);
      doglegSolver.set_max_radius(i);
      doglegSolver.set_initial_trust_radius(j);

      doglegSolver.solve(nonlinearProblem, solution);
      if(doglegSolver.has_converged())
        fileConvergence << doglegSolver.iterations()<< "\n";
      else
        fileConvergence << -1 << "\n";
    }
  }
  fileConvergence.close();

  system("python3 test_convergence.py");
  std::remove("convergence.dat");
  std::cout<<"-------------------------------"<<std::endl;
}


/**
 * @brief Test the convergence of the Newton method for different initial solutions. 
 *        This test method will only work for the input dimension 2.
 * 
 * @tparam NonlinearProblem 
 * @param nonlinearProblem 
 * @param domain 
 */
template<typename NonlinearProblem>
void testNewtonConvergence(const NonlinearProblem& nonlinearProblem, Vector<typename NonlinearProblem::number_type> domain){
  if(nonlinearProblem.size() != 2){
    HDNUM_ERROR("Error: This test method will only work for the input dimension 2");
  }
  typedef typename NonlinearProblem::number_type N;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"Test Newton method for different initial solutions"<< std::endl;
  std::cout<<"Domain: " << domain[0]  << " " << domain[1] << std::endl;
  Vector<N> solution(2);

  hdnum::Newton newton_solver;
  newton_solver.set_verbosity(0);
  newton_solver.set_linesearchsteps(50);
  newton_solver.set_maxit(1000);

  std::fstream fileConvergence("convergence.dat",std::ios::out);
  fileConvergence <<"1" << "\n"; 
  fileConvergence << domain[0] << " " << domain[1]<<"\n";

  //change the step size to quick up process
  N step_size = 0.2;
  for(N i=domain[0]; i<=domain[1]; i=i+step_size){
    for(N j=domain[0]; j<=domain[1]; j=j+step_size){
      solution[0] = i;
      solution[1] = j;

      newton_solver.solve(nonlinearProblem, solution);
      if(newton_solver.has_converged())
        fileConvergence << newton_solver.iterations()<< "\n";
      else
        fileConvergence << -1<< "\n";
    }
  }
  fileConvergence.close();

  system("python3 test_convergence.py");
  std::remove("convergence.dat");
  std::cout<<"-------------------------------"<<std::endl;

}

template<typename MinimizationProblem>
void testProjectedNewton(const MinimizationProblem& nonlinearProblem, DenseMatrix<typename MinimizationProblem::number_type> constraints, Vector<typename MinimizationProblem::number_type> lowerbound,  Vector<typename MinimizationProblem::number_type> upperbound){  
  typedef typename MinimizationProblem::number_type N;

   std::fstream file_loss("loss.dat",std::ios::out);

    N stepSize = 0.2;

    for(N i = -5; i<= 10 ; i = i+stepSize){
      for(N j = -5; j<= 10 ; j = j+stepSize){
        Vector<N> point = {i,j};
        Vector<N> result(1);
        nonlinearProblem.f(point, result);
        file_loss<< result[0]<< "\n";
      }
    }
    file_loss.close();
    
    std::fstream file_constraints("constraints.dat", std::ios::out);
    for(int i=0; i< constraints.rowsize(); ++i){
        file_constraints<<lowerbound[i]<<"   ";
        for(int j = 0; j< constraints.colsize();++j){
          file_constraints <<constraints(i,j) <<"   ";
        }
        file_constraints<<upperbound[i]<<"\n";
    }
    file_constraints.close();

  Vector<N> sol(2);
  sol[0] = 6;
  sol[1] = 2;
  ProjectedNewton proj;
  proj.solve(nonlinearProblem, sol, "a.dat");

  system("python3 test_projected_newton.py");

  std::remove("loss.dat");
  std::remove("constraints.dat");
  std::remove("a.dat");
}

int main ()
{
  // Examples of solving nonlinear models with the newton-raphson method:

  // Squared root problem: min_x x^2 - a. (Here, a=16)
  SquareRootProblem<Number> sqp(16);
  std::cout<< "Size of the model: " << sqp.size()<<std::endl;

  // Declare a Newton-Raphson solver and use default valued for the maximum nuber of iterations, 
  // the number of line search steps, the absolute limit for defect and the reduction factor.
  Newton newtonRapshon;
  // output the summary
  newtonRapshon.set_verbosity(1);

  // Solution of the nonlinear problem will be saved here.
  Vector<Number> x(1);

  x[0] = 25; //Change the initial value.
  
  newtonRapshon.solve(sqp,  x);

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  // Solution of the nonlinear problem will be saved here.
  Vector<Number> solution(2);

  // General nonlinear problems(For the loss function, you need to define a lambda function)
  auto problemRosenbrock = getNonlinearProblem(functionRosenbrock,solution);
  auto problemBohachesky = getNonlinearProblem(functionBohachesky, solution);
  auto problemBooth = getNonlinearProblem(functionBooth, solution);
  auto problemBranin = getNonlinearProblem(functionBranin, solution);
  auto problemMatyas = getNonlinearProblem(functionMatyas, solution);

  std::cout<< "Number of components of model: " << problemRosenbrock.size()<<std::endl;

  // Change settings of the solver
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
  ndlc.set_initial_trust_radius(100.0);
  ndlc.set_max_radius(100.0);
  ndlc.solve(problemRosenbrock,solution);
  std::cout<<"Solution: Newton Dog Cauchy: " << solution[0] << " " << solution[1] << std::endl;

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  ndlc.set_verbosity(0);
  newtonRapshon.set_verbosity(0);
  
  Vec sol(2);
  sol[0] = -5.0;
  sol[1] = -2.0;
  
  Vec domain = {-10.0,30.0};
  Vector<double> range = {0, 100};
  
  //testNewtonAgainstDogLeg(problemRosenbrock, sol, domain, 1, 1);
  //testDoglegConvergenceFixedRadius(problemBranin, 1.0,0.1, domain);
  //testDoglegConvergenceFixedInitialSolution(problemBranin,sol, range);
  //testNewtonConvergence(problemRosenbrock, domain);

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
  auto prob1= getNonlinearMinimizationProblem_Constrained(functionConstrained1, gradientConstrained1, constraints1, lowerbound1, upperbound1, s);

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
  auto prob2= getNonlinearMinimizationProblem_Constrained(functionConstrained2, gradientConstraiend2, constraints2, lowerbound2, upperbound2, s);

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
  auto prob3= hdnum::getNonlinearMinimizationProblem_Constrained(functionConstrained3, gradientConstraiend3, constraints3, lowerbound3, upperbound3, s);

  //proj.solve(prob3, s);

  testProjectedNewton(prob3, constraints3, lowerbound3, upperbound3);
}


