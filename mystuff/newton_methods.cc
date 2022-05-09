#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include "hdnum.hh"    // hdnum header
#include <vector>
#include <math.h>     
#include <map> 

#define PI 3.14159265

// Functions F(x), which describe the problem F(x) = 0

  auto FUNCTION_ROSENBROCK = [](hdnum::Vector<double> x) 
  { 
    hdnum::Vector<double> result(2);
    result[0] = -2* (1-x[0]) - 400 * (x[1] - x[0] * x[0]) * x[0];
    result[1] = 200 * (x[1] -x[0] * x[0]);

    return result;
  };

  auto FUNCTION_BOHACHEVSKY = [](hdnum::Vector<double> x) 
  { 
    hdnum::Vector<double> result(2);
    result[0] = 2 * x[0] + 0.3 * 3.0 * PI * sin(3.0 * PI * x[0]);
    result[1] = 4 * x[1] + 0.4 * 4.0 * PI * sin(4.0 * PI * x[1]);
    
    return result;
  };

  auto FUNCTION_BOOTH = [](hdnum::Vector<double> x) 
  { 
    hdnum::Vector<double> result(2);
    result[0] = 2 * (x[0] + 2 * x[1] - 7) + 4.0 * (2 * x[0] +x[1]-5);
    result[1] = 4 * (x[0] + 2 * x[1] - 7) + 2.0 * (2 * x[0] +x[1]-5);
    return result;
  };
  
  auto FUNCTION_BRANIN = [](hdnum::Vector<double> x) 
  { 
    hdnum::Vector<double> result(2);
    double a = 1;
    double b = 5.1 / ( 4.0 * PI * PI );
    double c = 5.0 / PI;

    double r = 6;
    double s = 10;
    double t = 1 / (8* PI);

    result[0] = 2 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r) * (-2 * b * x[0] + c) - s * (1 - t) * sin(x[0]);
    result[1] = 2 * a * (x[1] - b * x[0] * x[0] + c* x[0] -r);

    return result;
  };

/**
 * @brief Test Newton method against Newton dog leg cauchy method(with a fixed maximum trust radius and initial trust radius).
 *        This test method will only work for the input dimension 2.
 *        Both methods output a summary on the terminal and afterwards some results will be vizualized:
 *        - trajectory of the optimal solution(3D and 2D)
 *        - process of norm and reduction
 *        - vizualizing the trust radii of the Newton dog leg cauchy method
 * 
 * @tparam NonlinearProblem 
 * @param initial_solution 
 * @param nonlinearProblem 
 * @param domain 
 * @param max_trus_radius // default: 1.0
 * @param initial_trust_radius // default: 1.0
 */
template<typename NonlinearProblem>
void test_newton_against_dog_leg(hdnum::Vector<double> initial_solution, const NonlinearProblem& nonlinearProblem, hdnum::Vector<double> domain ,  double max_trust_radius= 1.0, double initial_trust_radius= 1.0){
    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<"Test Newton method against Newton Dog Leg Cauchy"<< std::endl;
    std::cout<<"Initial solution: " << initial_solution[0] << " " << initial_solution[1]<< std::endl;
    std::cout<<"Domain region: " << domain[0] << " " << domain[1] << std::endl;
    std::cout<<"Initial trust radius: " << initial_trust_radius << " , Maxmimum trust radius: " << max_trust_radius << std::endl;
    hdnum::Vector<double> solution = initial_solution;
    
    std::fstream file_loss("loss.dat",std::ios::out);
    file_loss << domain[0]<<" " << domain[1] << "\n";
    for(float i = domain[0]; i<= domain[1] ; i = i+0.2){
      for(float j = domain[0]; j<= domain[1] ; j = j+0.2){
        hdnum::Vector<double> point(2);
        point[0] = i;
        point[1] = j;
        hdnum::Vector<double> result(2);
        nonlinearProblem.F(point, result);
        double value = 0.5 * (result * result); 
        file_loss<< value << "\n";
      }
    }
    file_loss.close();

    hdnum::Newton newton_solver;
    newton_solver.set_verbosity(1);
    newton_solver.set_linesearchsteps(50);
    newton_solver.set_maxit(500);

    newton_solver.solve(nonlinearProblem, solution);
    newton_solver.save_data_in_file("newton_solver.dat");
    
    solution = initial_solution;

    hdnum::Newton_Dog_Leg_Cauchy dog_leg_solver;
    dog_leg_solver.set_verbosity(1);
    dog_leg_solver.set_maxit(500);
    dog_leg_solver.set_reduction(1e-15);
    dog_leg_solver.set_initial_trust_radius(initial_trust_radius);
    dog_leg_solver.set_max_radius(max_trust_radius);

    dog_leg_solver.solve(nonlinearProblem, solution);
    dog_leg_solver.save_data_in_file("dog_leg_solver.dat");

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
 * @param initial_trust_radius 
 * @param max_trust_radius 
 * @param nonlinearProblem 
 * @param domain 
 */
template<typename NonlinearProblem>
void test_dog_leg_convergence_fixed_radius(double initial_trust_radius, double max_trust_radius, const NonlinearProblem& nonlinearProblem, hdnum::Vector<double> domain){
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"Test Newton dog leg cauchy method for different initial solutions and a fixed initial and maximum trust radius"<< std::endl;
  std::cout<<"Initial trust radius: " << initial_trust_radius << " , maximum trust radius: " << max_trust_radius << std::endl;
  std::cout<<"Domain: " << domain[0] << " " << domain[1] << std::endl;

  hdnum::Vector<double> solution(nonlinearProblem.size());
  hdnum::Newton_Dog_Leg_Cauchy dog_leg_solver;
  dog_leg_solver.set_verbosity(0);
  dog_leg_solver.set_maxit(1000);
  dog_leg_solver.set_max_radius(max_trust_radius);
  dog_leg_solver.set_initial_trust_radius(initial_trust_radius);

  std::fstream file_convergence("convergence.dat",std::ios::out);
  file_convergence <<"1" << "\n"; 
  file_convergence << domain[0] << " " << domain[1]<<"\n";

  //change the step size to quick up process
  double step_size = 0.2;
  double value = (domain[1] - domain[0]) / step_size;

  for(float i=domain[0]; i<=domain[1]; i=i+step_size){
    for(float j=domain[0]; j<=domain[1]; j=j+step_size){
      solution[0] = i;
      solution[1] = j;

      dog_leg_solver.solve(nonlinearProblem, solution);
      if(dog_leg_solver.has_converged())
        file_convergence << dog_leg_solver.iterations()<< "\n";
      else
        file_convergence << -1<< "\n";
    }
  }
  
  file_convergence.close();

  system("python3 test_convergence.py");
  std::remove("convergence.dat");
  std::cout<<"-------------------------------"<<std::endl;
}

/**
 * @brief Test the convergence of the Newton dog leg cauchy method for different initial and maximum trust radii 
 *        and a fixed initial solution. This method will work for any input dimension.
 * 
 * @tparam NonlinearProblem 
 * @param initial_solution 
 * @param nonlinearProblem 
 * @param range // [minimum max_trust_radius, maximum max_trust_radius]
 */
template<typename NonlinearProblem>
void test_dog_leg_convergence_fixed_initial_solution(hdnum::Vector<double> initial_solution, const NonlinearProblem& nonlinearProblem, hdnum::Vector<double> range){
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"Test Newton dog leg cauchy method for different initial and maximum trust radii and a fixed initial solution"<< std::endl;
  std::cout<<"Initial solution: ";
  for(auto sol: initial_solution){
    std::cout<< sol;
  }
  std::cout<<std::endl;
  std::cout<<"Range: " << range[0] << " " << range[1] << std::endl;

  hdnum::Vector<double> solution(2);
  hdnum::Newton_Dog_Leg_Cauchy dog_leg_solver;

  std::fstream file_convergence("convergence.dat",std::ios::out);
  file_convergence <<"2" << "\n"; 
  file_convergence << range[0] << " " << range[1]<<"\n";

  //change the step size to quick up process
  double step_size = 0.2;

  for(float i=range[0]; i<=range[1]; i=i+step_size){
    for(float j=range[0]; j<=range[1]; j=j+step_size){
      if(j > i){
        file_convergence << -1<< "\n";
        continue;
      }
      solution = initial_solution;

      dog_leg_solver.set_verbosity(0);
      dog_leg_solver.set_maxit(1000);
      dog_leg_solver.set_max_radius(i);
      dog_leg_solver.set_initial_trust_radius(j);

      dog_leg_solver.solve(nonlinearProblem, solution);
      if(dog_leg_solver.has_converged())
        file_convergence << dog_leg_solver.iterations()<< "\n";
      else
        file_convergence << -1 << "\n";
    }
  }
  file_convergence.close();

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
void test_newton_convergence(const NonlinearProblem& nonlinearProblem, hdnum::Vector<double> domain){
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"Test Newton method for different initial solutions"<< std::endl;
  std::cout<<"Domain: " << domain[0]  << " " << domain[1] << std::endl;
  hdnum::Vector<double> solution(2);

  hdnum::Newton newton_solver;
  newton_solver.set_verbosity(0);
  newton_solver.set_linesearchsteps(50);
  newton_solver.set_maxit(1000);

  std::fstream file_convergence("convergence.dat",std::ios::out);
  file_convergence <<"1" << "\n"; 
  file_convergence << domain[0] << " " << domain[1]<<"\n";

  //change the step size to quick up process
  double step_size = 0.2;
  for(float i=domain[0]; i<=domain[1]; i=i+step_size){
    for(float j=domain[0]; j<=domain[1]; j=j+step_size){
      solution[0] = i;
      solution[1] = j;

      newton_solver.solve(nonlinearProblem, solution);
      if(newton_solver.has_converged())
        file_convergence << newton_solver.iterations()<< "\n";
      else
        file_convergence << -1<< "\n";
    }
  }
  file_convergence.close();

  system("python3 test_convergence.py");
  std::remove("convergence.dat");
  std::cout<<"-------------------------------"<<std::endl;

}



int main ()
{
  // Examples of solving nonlinear models with the newton-raphson method:

  // Squared root problem: min_x x^2 - a. (Here, a=16)
  hdnum::SquareRootProblem<double> sqp(16);
  std::cout<< "Size of the model: " << sqp.size()<<std::endl;

  // Declare a Newton-Raphson solver and use default valued for the maximum nuber of iterations, 
  // the number of line search steps, the absolute limit for defect and the reduction factor.
  hdnum::Newton newton_raphson;
  // output the summary
  newton_raphson.set_verbosity(1);

  // Solution of the nonlinear problem will be saved here.
  hdnum::Vector<double> x(1);

  x[0] = 25; //Change the initial value.
  
  newton_raphson.solve(sqp,  x);

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  // Solution of the nonlinear problem will be saved here.
  hdnum::Vector<double> solution(2);

  // General nonlinear problems(For the loss function, you need to define a lambda function)
  auto PROBLEM_ROSENBROCK = hdnum::getNonlinearProblem(FUNCTION_ROSENBROCK,solution);
  auto PROBLEM_BOHACHEVSKY = hdnum::getNonlinearProblem(FUNCTION_BOHACHEVSKY, solution);
  auto PROBLEM_BOOTH = hdnum::getNonlinearProblem(FUNCTION_BOOTH, solution);
  auto PROBLEM_BRANIN = hdnum::getNonlinearProblem(FUNCTION_BRANIN, solution);

  std::cout<< "Number of components of model: " << PROBLEM_ROSENBROCK.size()<<std::endl;

  // Change settings of the solver
  newton_raphson.set_linesearchsteps(50);
  newton_raphson.set_maxit(10000);
  newton_raphson.set_reduction(1e-15);

  solution[0] = -4;
  solution[1] = -4;

  newton_raphson.solve(PROBLEM_ROSENBROCK, solution);
  std::cout<< "Solution Newton: " << solution[0] << " " << solution[1] << std::endl;

  //newton_raphson.save_data_in_file("problem_rosenbrock_newton_initial_point_2_2.dat");

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  hdnum::Newton_Dog_Leg_Cauchy ntr;
  solution[0] = -4.0;
  solution[1] = -4.0;
  ntr.set_verbosity(1);
  ntr.set_maxit(1000);
  ntr.set_reduction(1e-15);
  ntr.set_initial_trust_radius(1.0);
  ntr.set_max_radius(1.0);
  ntr.solve(PROBLEM_ROSENBROCK,solution);
  std::cout<<"Solution: Newton Dog Cauchy: " << solution[0] << " " << solution[1] << std::endl;

  std::cout<<"-----------------------------------------------------------"<<std::endl;

  //ntr.save_data_in_file("problem_rosenbrock_dog_leg_initial_point_2_2.dat");

  ntr.set_verbosity(0);
  newton_raphson.set_verbosity(0);

  
  hdnum::Vector<double> sol(2);
  sol[0] = 5.0;
  sol[1] = 10.0;
  
  hdnum::Vector<double> domain = {-6.0,15.0};
  hdnum::Vector<double> range = {1.0, 100};
  
  //test_newton_against_dog_leg(sol, PROBLEM_BRANIN, domain, 5.0, 5.0);
  //test_dog_leg_convergence_fixed_radius(1.0,0.1, PROBLEM_BRANIN, domain);
  //test_dog_leg_convergence_fixed_initial_solution(sol, PROBLEM_BRANIN, range);
  //test_newton_convergence(PROBLEM_BRANIN, domain);
}


