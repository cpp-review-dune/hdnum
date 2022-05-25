#ifndef NEWTONVISUALIZATION
#define NEWTONVISUALIZATION

#include <iostream>
#include "hdnum.hh"

using namespace hdnum;

std::string pythonScriptProjectedNewton = "testProjectedNewton.py";
std::string pythonScriptNewtonAgainstDogleg = "testNewtonAgainstDogLeg.py";
std::string pythonScriptConvergence = "testConvergence.py";

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
    std::cout<<"Domain: " << domain[0] << " " << domain[1] << std::endl;
    std::cout<<"Initial trust radius: " << initialTrustRadius << " , Maxmimum trust radius: " << maxTrustRadius << std::endl;
    Vector<N> solution = initialSolution;
    
    std::fstream file_loss("loss.dat",std::ios::out); // create file to save loss(x) = 0.5 * F(x)^T * F(x)
    file_loss << domain[0]<<" " << domain[1] << "\n";

    N stepSize = 0.2;

    // compute losses over the domain
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
    doglegSolver.set_maxit(500); // set maximum number of iteration
    doglegSolver.set_reduction(1e-15);
    doglegSolver.setInitialTrustRadius(initialTrustRadius);
    doglegSolver.setMaxTrustRadius(maxTrustRadius);

    doglegSolver.solve(nonlinearProblem, solution,"dog_leg_solver.dat");

    // call python script to visulize results
    std::string command = "python3 " + pythonScriptNewtonAgainstDogleg; 
    system(command.c_str());

    // delete files
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

  NewtonDogLegCauchy doglegSolver;
  doglegSolver.set_verbosity(0);
  doglegSolver.set_maxit(1000);
  doglegSolver.setMaxTrustRadius(maxTrustRadius);
  doglegSolver.setInitialTrustRadius(initialTrustRadius);

  std::fstream fileConvergence("convergence.dat",std::ios::out); // saves results of type: has converged or not
  fileConvergence <<"1" << "\n"; 
  fileConvergence << domain[0] << " " << domain[1]<<"\n";

  std::fstream fileConvergence2("convergence2.dat",std::ios::out); // saves results of type: to which optimum did the method converge
  std::vector<std::vector<N>> solutions; // set of solutions
  Vector<N> solution(nonlinearProblem.size());

  N step_size = 0.2;  // change the step size to quick up process

  // check convergence over the domain
  for(N i=domain[0]; i<=domain[1]; i=i+step_size){
    for(N j=domain[0]; j<=domain[1]; j=j+step_size){
      solution[0] = i;
      solution[1] = j;

      doglegSolver.solve(nonlinearProblem, solution);
      std::vector<N> solution_temp = {round(solution[0]), round(solution[1])};
      if(doglegSolver.has_converged()){
        fileConvergence << doglegSolver.iterations()<< "\n";

        bool found = false;
        for(int k = 0; k< solutions.size(); ++k){
          if(solutions[k][0] == solution_temp[0] && solutions[k][1] == solution_temp[1]){ // checks if solution is already in the set
            fileConvergence2<<k<<"\n"; // save the set index
            found = true;
            break;
          }
        }
        if(!found){
          solutions.push_back(solution_temp); // push solution to the set
          fileConvergence2<<solutions.size() -1<<"\n";
        }
      }
      else{
        fileConvergence << -1 << "\n";
        fileConvergence2 << -1 << "\n";
      }
    }
  }
  
  fileConvergence.close();
  fileConvergence2.close();

  // call python script to visulize results
  std::string command = "python3 " + pythonScriptConvergence; 
  system(command.c_str());

  // delete file
  std::remove("convergence.dat");
  std::remove("convergence2.dat");
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

  NewtonDogLegCauchy doglegSolver;

  std::fstream fileConvergence("convergence.dat",std::ios::out); // saves results of type: has converged or not
  fileConvergence <<"2" << "\n"; 
  fileConvergence << range[0] << " " << range[1]<<"\n";

  std::fstream fileConvergence2("convergence2.dat",std::ios::out); // saves results of type: to which optimum did the method converge
  std::vector<std::vector<N>> solutions; // set of solutions
  Vector<N> solution(nonlinearProblem.size());

  double stepSize = 0.2;  // change the step size to quick up process

  // check convergence over the domain
  for(double i=range[0]; i<=range[1]; i=i+stepSize){
    for(double j=range[0]; j<=range[1]; j=j+stepSize){
      if(j > i){
        fileConvergence << -1 << "\n";
        fileConvergence2 << -1 << "\n";
        continue;
      }
      solution = initialSolution;

      doglegSolver.set_verbosity(0);
      doglegSolver.set_maxit(500);
      doglegSolver.setMaxTrustRadius(i);
      doglegSolver.setInitialTrustRadius(j);

      doglegSolver.solve(nonlinearProblem, solution);
      std::vector<N> solution_temp = {round(solution[0]), round(solution[1])};
      if(doglegSolver.has_converged()){
        fileConvergence << doglegSolver.iterations()<< "\n";
        bool found = false;
        for(int k = 0; k< solutions.size(); ++k){
          if(solutions[k][0] == solution_temp[0] && solutions[k][1] == solution_temp[1]){ // checks if solution is already in the set
            fileConvergence2<<k<<"\n"; // save the set index
            found = true;
            break;
          }
        }
        if(!found){
          solutions.push_back(solution_temp); // push solution to the set
          fileConvergence2<<solutions.size() -1<<"\n";
        }
      }
      else{
        fileConvergence << -1 << "\n";
        fileConvergence2 << -1 << "\n";
      }
    }
  }
  fileConvergence.close();
  fileConvergence2.close();

  // call python script to visulize results
  std::string command = "python3 " + pythonScriptConvergence;
  system(command.c_str());

  // delete file
  std::remove("convergence.dat");
  std::remove("convergence2.dat");
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

  hdnum::Newton newton_solver;
  newton_solver.set_verbosity(0);
  newton_solver.set_linesearchsteps(50);
  newton_solver.set_maxit(1000);

  std::fstream fileConvergence("convergence.dat",std::ios::out); // saves results of type: has converged or not
  fileConvergence <<"1" << "\n"; 
  fileConvergence << domain[0] << " " << domain[1]<<"\n";

  std::fstream fileConvergence2("convergence2.dat",std::ios::out); // saves results of type: to which optimum did the method converge

  N step_size = 0.2; // change the step size to quick up process

  std::vector<std::vector<N>> solutions; // set of solutions
  Vector<N> solution(nonlinearProblem.size());

  // check convergence over the domain
  for(N i=domain[0]; i<=domain[1]; i=i+step_size){
    for(N j=domain[0]; j<=domain[1]; j=j+step_size){
      solution[0] = i;
      solution[1] = j;

      newton_solver.solve(nonlinearProblem, solution);
      std::vector<N> solution_temp = {round(solution[0]), round(solution[1])};
      if(newton_solver.has_converged()){
        fileConvergence << newton_solver.iterations()<< "\n";
        bool found = false;
        for(int k = 0; k< solutions.size(); ++k){
          if(solutions[k][0] == solution_temp[0] && solutions[k][1] == solution_temp[1]){ // checks if solution is already in the set
            fileConvergence2<<k<<"\n"; // save the set index
            found = true;
            break;
          }
        }
        if(!found){
          solutions.push_back(solution_temp); // push solution to the set
          fileConvergence2<<solutions.size() -1<<"\n";
        }
      }
      else{
        fileConvergence << -1 << "\n";
        fileConvergence2 << -1 << "\n";
      }
    }
  }
  fileConvergence.close();
  fileConvergence2.close();

  // call python script to visulize results
  std::string command = "python3 " + pythonScriptConvergence;
  system(command.c_str());

  // delete file
  std::remove("convergence.dat");
  std::remove("convergence2.dat");
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

  std::string command = "python3 " + pythonScriptProjectedNewton;
  system(command.c_str());

  std::remove("loss.dat");
  std::remove("constraints.dat");
  std::remove("a.dat");
}

#endif