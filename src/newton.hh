// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_NEWTON_HH
#define HDNUM_NEWTON_HH

#include "lr.hh"
#include <type_traits>
#include <vector>
#include <stdlib.h>     /* srand, rand */

/** @file
 *  @brief Newton's method with line search
 */

namespace hdnum {

  /** @brief Example class for a nonlinear model F(x) = 0;

      This example solves F(x) = x*x - a = 0

      \tparam N a type representing x and F components
  */
  template<class N>
  class SquareRootProblem
  {
  public:
    /** \brief export size_type */
    typedef std::size_t size_type;

    /** \brief export number_type */
    typedef N number_type;

    //! constructor stores parameter lambda
    SquareRootProblem (number_type a_)
      : a(a_)
    {}

    //! return number of componentes for the model
    std::size_t size () const
    {
      return 1;
    }

    //! model evaluation
    void F (const Vector<N>& x, Vector<N>& result) const
    {
      result[0] = x[0]*x[0] - a;
    }

    //! jacobian evaluation needed for implicit solvers
    void F_x (const Vector<N>& x, DenseMatrix<N>& result) const
    {
      result[0][0] = number_type(2.0)*x[0];
    }

  private:
    number_type a;
  };


  /** @brief A generic problem class that can be set up with a lambda defining F(x)=0

      \tparam Lambda mapping a Vector to a Vector
      \tparam Vec    the type for the Vector
  */
  template<typename Lambda, typename Vec>
  class GenericNonlinearProblem
  {
    Lambda lambda; // lambda defining the problem "lambda(x)=0"
    size_t s;
    typename Vec::value_type eps;
  
  public:
    /** \brief export size_type */
    typedef std::size_t size_type;

    /** \brief export number_type */
    typedef typename Vec::value_type number_type;

    //! constructor stores parameter lambda
    GenericNonlinearProblem (const Lambda& l_, const Vec& x_, number_type eps_ = 1e-7)
      : lambda(l_), s(x_.size()), eps(eps_)
    {}

    //! return number of componentes for the model
    std::size_t size () const
    {
      return s;
    }

    //! model evaluation
    void F (const Vec& x, Vec& result) const
    {
      result = lambda(x);
    }

    //! jacobian evaluation needed for implicit solvers
    void F_x (const Vec& x, DenseMatrix<number_type>& result) const
    {
      Vec Fx(x.size());
      F(x,Fx);
      Vec z(x);
      Vec Fz(x.size());
    
      // numerische Jacobimatrix
      for (int j=0; j<result.colsize(); ++j)
        {
          auto zj = z[j];
          auto dz = (1.0+abs(zj))*eps;
          z[j] += dz;
          F(z,Fz);
          for (int i=0; i<result.rowsize(); i++)
            result[i][j] = (Fz[i]-Fx[i])/dz;
          z[j] = zj;
        }
    }
  };

  /** @brief A function returning a problem class

      Automatic template parameter extraction makes fiddling with types unnecessary.

      \tparam F  a lambda mapping a Vector to a Vector
      \tparam X  the type for the Vector
  */
  template<typename F, typename X>
  GenericNonlinearProblem<F,X> getNonlinearProblem (const F& f, const X& x, typename X::value_type eps = 1e-7)
  {
    return GenericNonlinearProblem<F,X>(f,x,eps);
  }


  /**
   * @brief A generic minimization problem class that can be set up with a lambda(Objective) definining min F(x),
   *        a lambda(gradient) defining the gradient of F and constraints defined by lower/upper bounds and a transformation matrix A.
   * 
   * @tparam Objective F(x) 
   * @tparam Gradient gradient of F
   * @tparam Vec the type of the Vector
   */
  template<typename Objective, typename Gradient, typename Vec>
  class GenericNonlinearMinimizationProblem_Constrained
  {
    Objective objective; 
    Gradient gradient;

    size_t s;
    typename Vec::value_type eps;
  
  public:
    /** \brief export size_type */
    typedef std::size_t size_type;

    /** \brief export number_type */
    typedef typename Vec::value_type number_type;

    DenseMatrix<number_type> A;
    Vector<number_type> upperBounds;
    Vector<number_type> lowerBounds;

    GenericNonlinearMinimizationProblem_Constrained (const Objective& o, const Gradient& g, const DenseMatrix<double> constraints, const Vector<double> lb, const Vector<double> up, const  Vec& x_, number_type eps_ = 1e-7)
      : objective(o), gradient(g), s(x_.size()), eps(eps_), upperBounds(up), lowerBounds(lb), A(constraints)
    {}

    //! return number of componentes for the model
    std::size_t size () const
    {
      return s;
    }

    // objective
    void f(const Vec& x, Vec& result) const
    {
      result = objective(x);
    }

    // gradient
    void g(const Vec& x, Vec& result) const
    {
      result = gradient(x);
    }

    // Hessian
    void H(const Vec& x, DenseMatrix<number_type>& result) const
    {
      Vec Fx(x.size());
      g(x,Fx);
      Vec z(x);
      Vec Fz(x.size());
    
      for (int j=0; j<result.colsize(); ++j)
        {
          auto zj = z[j];
          auto dz = (1.0+abs(zj))*eps;
          z[j] += dz;
          g(z,Fz);
          for (int i=0; i<result.rowsize(); i++)
            result[i][j] = (Fz[i]-Fx[i])/dz;
          z[j] = zj;
        }
    }
    
    void project_to_feasible_space(Vec& x) const{
      for(int i = 0; i<x.size(); ++i){
        if(x[i] < lowerBounds[i]){
          x[i] = lowerBounds[i];
        }
        else if(x[i] > upperBounds[i]){
          x[i] = upperBounds[i];
        }
      }
    }

    int constraintViolated(Vec& x) const{
      auto y = A * x;
      for(int i=0; i<y.size();++i){
        if(y[i] < lowerBounds[i] || y[i] > upperBounds[i]){
          return i;
        }
      }

      return -1;
    }
  };

  /**
   * @brief A function returning a minimization problem class with constraints
   * 
   * @tparam Objective F(x) 
   * @tparam Gradient gradient of F
   * @tparam X the type of the Vector
   * @param constraints transformation matrix A
   * @param lb lower bound
   * @param up upper bound
   * @return GenericNonlinearMinimizationProblem_Constrained<Objective,Gradient,X> 
   */
  template<typename Objective,typename Gradient, typename X>
  GenericNonlinearMinimizationProblem_Constrained<Objective,Gradient,X> getNonlinearMinimizationProblem_Constrained(const Objective& o, const Gradient& g, const DenseMatrix<double> constraints, const Vector<double> lb, const Vector<double> up, const  X& x, typename X::value_type eps = 1e-7)
  {
    return GenericNonlinearMinimizationProblem_Constrained<Objective, Gradient,X>(o,g,constraints, lb, up, x,eps);
  }

  /** @brief Solve nonlinear problem using a damped Newton method

      The Newton solver is parametrized by a model. The model also
      exports all relevant types for types.

  */
  class Newton
  {

  private: 
    typedef std::size_t size_type;

    template<typename N, typename Real>
    void saveDataInFile(std::string filename, std::vector<std::vector<N>> iterationPoints, std::vector<std::vector<N>> directions, std::vector<Real> stepSizes, std::vector<Real> norms, std::vector<Real> reductions,std::vector<Real> losses){
        std::fstream file(filename.c_str(),std::ios::out);
        file<< "#This file contains the outputs from the newton-raphson method solver. The outputs are ordered in the following way: \n";
        file<< "#Iteration point (X,Y)) | directions | step size | norm | reduction \n";
        for(int i = 0; i< iterations_taken;++i){
            Vector<N> x(2);
            x[0] = iterationPoints[i][0]; 
            x[1] = iterationPoints[i][1];

            Vector<N> direction(2);
            direction[0] = directions[i][0];
            direction[1] = directions[i][1];
            file << x[0] << "   " << x[1] << "   " << direction[0]<< "   " << direction[1] << "   "<< stepSizes[i] << "   " << norms[i] << "   " << reductions[i]<< "   " << losses[i] << "\n"; 

        }
    }

  public:
    //! constructor stores reference to the model
    Newton ()
      : maxit(25), linesearchsteps(10), verbosity(0), 
        reduction(1e-14), abslimit(1e-30), converged(false)
    {}

    //! maximum number of iterations before giving up
    void set_maxit (size_type n)
    {
      maxit = n;
    }

    void set_sigma (double sigma_)
    {
    
    }

    //! maximum number of steps in linesearch before giving up
    void set_linesearchsteps (size_type n)
    {
      linesearchsteps = n;
    }

    //! control output given 0=nothing, 1=summary, 2=every step, 3=include line search
    void set_verbosity (size_type n)
    {
      verbosity = n;
    }

    //! basolute limit for defect
    void set_abslimit (double l)
    {
      abslimit = l;
    }

    //! reduction factor
    void set_reduction (double l)
    {
      reduction = l;
    }


    //! do one step
    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x, std::string filename = "") 
    {
      //clear_data();

      typedef typename M::number_type N;
      // In complex case, we still need to use real valued numbers for residual norms etc.
      using Real = typename std::conditional<std::is_same<std::complex<double>, N>::value, double, N>::type;

      std::vector<std::vector<N>> iteration_points;
      std::vector<std::vector<N>> directions;
      std::vector<Real> step_sizes;
      std::vector<Real> norms;
      std::vector<Real> reductions;
      std::vector<Real> loss;


      Vector<N> r(model.size());              // residual
      DenseMatrix<N> A(model.size(),model.size()); // Jacobian matrix
      Vector<N> y(model.size());              // temporary solution in line search
      Vector<N> z(model.size());              // solution of linear system
      Vector<N> s(model.size());              // scaling factors
      Vector<size_type> p(model.size());                 // row permutations
      Vector<size_type> q(model.size());                 // column permutations

      model.F(x,r);                                     // compute nonlinear residualz

      Real R0(std::abs(norm(r)));                          // norm of initial residual
      Real R(R0);                                // current residual norm
      if (verbosity>=1)
        {
          std::cout << "Newton " 
                    << "   norm=" << std::scientific << std::showpoint 
                    << std::setprecision(4) << R0
                    << std::endl;
        }

      converged = false;
      for (size_type i=1; i<=maxit; i++)                // do Newton iterations
        {
          iteration_points.push_back( {x[0], x[1]});
          norms.push_back(std::abs(norm(r)));
          loss.push_back(0.5 * std::abs(norm(r)) * std::abs(norm(r)));

          // check absolute size of residual
          if (R<=abslimit)
            {
              converged = true;
              if(filename.size()>0){
                saveDataInFile(filename,iteration_points,directions, step_sizes, norms, reductions, loss);
              }
              return;
            } 

          // solve Jacobian system for update
          model.F_x(x,A);                               // compute Jacobian matrix

          row_equilibrate(A,s);                         // equilibrate rows

          lr_fullpivot(A,p,q);                          // LR decomposition of A

          z = N(0.0);                                   // clear solution
          apply_equilibrate(s,r);                       // equilibration of right hand side

          permute_forward(p,r);                         // permutation of right hand side
          solveL(A,r,r);                                // forward substitution
          solveR(A,z,r);                                // backward substitution
          permute_backward(q,z);                        // backward permutation

          directions.push_back({-z[0], -z[1]});

          // line search
          Real lambda(1.0);                      // start with lambda=1
          for (size_type k=0; k<linesearchsteps; k++)
            {
              y = x;                              
              y.update(-lambda,z);                       // y = x+lambda*z
              model.F(y,r);                             // r = F(y)
              Real newR(std::abs(norm(r)));                // compute norm
              if (verbosity>=3)
                {
                  std::cout << "    line search "  << std::setw(2) << k 
                            << " lambda=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << lambda
                            << " norm=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << newR
                            << " red=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << newR/R
                            << std::endl;
                }
              if (newR<(1.0-0.25*lambda)*R)            // check convergence
                {
                  step_sizes.push_back(lambda);
                  reductions.push_back(newR/R);
                  if (verbosity>=2)
                    {
                      std::cout << "  step"  << std::setw(3) << i 
                                << " norm=" << std::scientific << std::showpoint 
                                << std::setprecision(4) << newR
                                << " red=" << std::scientific << std::showpoint 
                                << std::setprecision(4) << newR/R
                                << std::endl;
                    }
                  x = y;
                  R = newR;
                  break;                                // continue with Newton loop
                }
              else lambda *= 0.5;                       // reduce damping factor
              if (k==linesearchsteps-1)
                {
                  if (verbosity>=3){
                    std::cout << "    line search not converged within " << linesearchsteps << " steps" << std::endl; 
                    if(filename.size()>0){
                        saveDataInFile(filename,iteration_points,directions, step_sizes, norms, reductions, loss);
                    }
              }

                  return;
                }
            }

          // check convergence
          if (R<=reduction*R0)
            {
              if (verbosity>=1)
                {
                  std::cout << "Newton converged in "  << i << " steps"
                            << " reduction=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << R/R0
                            << std::endl;
                }
              iterations_taken = i;
              converged = true;
              if(filename.size()>0){
                saveDataInFile(filename,iteration_points,directions, step_sizes, norms, reductions, loss);
              }
              return;
            }
          if (i==maxit)
            {
              iterations_taken = i;
              if (verbosity>=1){
                std::cout << "Newton not converged within " << maxit << " iterations" << std::endl;
                                    if(filename.size()>0){
                        saveDataInFile(filename,iteration_points,directions, step_sizes, norms, reductions, loss);
                    }
              }
            }
        }
    }


    
    bool has_converged () const
    {
      return converged;
    }
    size_type iterations() const {
      return iterations_taken;
    }

  private:
    size_type maxit;
    mutable size_type iterations_taken = -1;
    size_type linesearchsteps;
    size_type verbosity;
    double reduction;
    double abslimit;
    mutable bool converged;
  };

 



  /** @brief Solve nonlinear problem using a fixed point iteration

      solve F(x) = 0.

      \f[ x = x - \sigma*F(x) \f]

  */
  class Banach
  {
    typedef std::size_t size_type;

  public:
    //! constructor stores reference to the model
    Banach ()
      : maxit(25), linesearchsteps(10), verbosity(0), 
        reduction(1e-14), abslimit(1e-30),  sigma(1.0), converged(false)
    {}

    //! maximum number of iterations before giving up
    void set_maxit (size_type n)
    {
      maxit = n;
    }

    //! damping parameter
    void set_sigma (double sigma_)
    {
      sigma = sigma_;
    }

    //! maximum number of steps in linesearch before giving up
    void set_linesearchsteps (size_type n)
    {
      linesearchsteps = n;
    }

    //! control output given 0=nothing, 1=summary, 2=every step, 3=include line search
    void set_verbosity (size_type n)
    {
      verbosity = n;
    }

    //! basolute limit for defect
    void set_abslimit (double l)
    {
      abslimit = l;
    }

    //! reduction factor
    void set_reduction (double l)
    {
      reduction = l;
    }

    //! do one step
    template<class M>
    void solve (const M& model, Vector<typename M::number_type>& x) const
    {
      typedef typename M::number_type N;
      Vector<N> r(model.size());              // residual
      Vector<N> y(model.size());              // temporary solution in line search

      model.F(x,r);                           // compute nonlinear residual
      N R0(norm(r));                          // norm of initial residual
      N R(R0);                                // current residual norm
      if (verbosity>=1)
        {
          std::cout << "Banach " 
                    << " norm=" << std::scientific << std::showpoint 
                    << std::setprecision(4) << R0
                    << std::endl;
        }

      converged = false;
      for (size_type i=1; i<=maxit; i++)                // do iterations
        {
          // check absolute size of residual
          if (R<=abslimit)
            {
              converged = true;
              return;
            } 

          // next iterate
          y = x;                                    
          y.update(-sigma,r);                       // y = x+lambda*z
          model.F(y,r);                             // r = F(y)
          N newR(norm(r));                // compute norm
          if (verbosity>=2)
            {
              std::cout << "    "  << std::setw(3) << i 
                        << " norm=" << std::scientific << std::showpoint 
                        << std::setprecision(4) << newR
                        << " red=" << std::scientific << std::showpoint 
                        << std::setprecision(4) << newR/R
                        << std::endl;
            }
          x = y;                                // accept new iterate
          R = newR;                             // remember new norm

          // check convergence
          if (R<=reduction*R0 || R<=abslimit)
            {
              if (verbosity>=1)
                {
                  std::cout << "Banach converged in "  << i << " steps"
                            << " reduction=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << R/R0
                            << std::endl;
                }
              converged = true;
              return;
            }
        }
    }
    
    bool has_converged () const
    {
      return converged;
    }

  private:
    size_type maxit;
    size_type linesearchsteps;
    size_type verbosity;
    double reduction;
    double abslimit;
    double sigma;
    mutable bool converged;
  };

/** @brief  Solve nonlinear problem F(x) = 0 by defining an eqivalent minimizing problem arg min_x f(x) = 0.5 * F(x)^T * F(x).
 * (Uses Trust region method where a subproblem is solved by different direction choices:
 *  Newton direction
 *  Steepest descent direction
 *  Dog leg direction)
 */
  class NewtonDogLegCauchy{
    private:
    
    typedef std::size_t size_type;

    template<typename N, typename Real>
    void saveDatainFile(std::string filename, std::vector<std::vector<N>> iteration_points, std::vector<std::vector<N>> newton_directions, std::vector<std::vector<N>> steepest_descent_directions, std::vector<std::vector<N>> dog_leg_directions, std::vector<double> trust_radius_list, std::vector<Real> norms, std::vector<Real> reductions, std::vector<Real> loss){
      std::fstream file(filename.c_str(),std::ios::out);
      file<< "#This file contains the outputs from the newton-dog-leg-cauchy method solver. The outputs are ordered in the following way: \n";
      file<< "#Iteration point (X,Y)) | directions(Newton-Steepest-Dog leg) | trust radius | norm | reduction \n";
      for(int i = 0; i< iterations_taken;++i){
        Vector<N> x(2);
        x[0] = iteration_points[i][0]; 
        x[1] = iteration_points[i][1];

        Vector<N> newton_direction(2);
        newton_direction[0] = newton_directions[i][0];
        newton_direction[1] = newton_directions[i][1];

        Vector<N> steepest_descent_direction(2);
        steepest_descent_direction[0] = steepest_descent_directions[i][0];
        steepest_descent_direction[1] = steepest_descent_directions[i][1];

        Vector<N> dog_leg_direction(2);
        dog_leg_direction[0] = dog_leg_directions[i][0];
        dog_leg_direction[1] = dog_leg_directions[i][1];

        file << x[0] << "   " << x[1] << "   ";
        file << newton_direction[0] << "   " << newton_direction[1] << "   ";
        file << steepest_descent_direction[0] << "   " << steepest_descent_direction[1] << "   ";
        file << dog_leg_direction[0] << "   " << dog_leg_direction[1] << "   ";
        file << trust_radius_list[i] << "   ";
        file << norms[i] << "   ";
        file << reductions[i] << "   ";
        file << loss[i] <<"\n";
      }
    }

    public:
      
    NewtonDogLegCauchy ()
      : maxit(25), verbosity(0), 
        reduction(1e-14), abslimit(1e-30), converged(false)
    {}

    void set_maxit (size_type n)
    {
      maxit = n;
    }

    //! basolute limit for defect
    void set_abslimit (double l)
    {
      abslimit = l;
    }

    //! control output given 0=nothing, 1=summary, 2=every step 
    void set_verbosity (size_type n)
    {
      verbosity = n;
    }

    //! reduction factor
    void set_reduction (double l)
    {
      reduction = l;
    }

    void set_max_radius(double max){
      max_radius = max;
    }

    void set_initial_trust_radius(double rad){
      initial_trust_radius = rad;
    }

    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x, std::string filename="") 
    {
      //clear_data();

      typedef typename M::number_type N;
      // In complex case, we still need to use real valued numbers for residual norms etc.
      using Real = typename std::conditional<std::is_same<std::complex<double>, N>::value, double, N>::type;

      std::vector<std::vector<N>> iteration_points;
      std::vector<std::vector<N>> newton_directions;
      std::vector<std::vector<N>> steepest_descent_directions;
      std::vector<std::vector<N>> dog_leg_directions;

      std::vector<double> trust_radius_list;
      std::vector<Real> norms;
      std::vector<Real> reductions;
      std::vector<Real> loss;

      Vector<N> F(model.size());              // residual
      DenseMatrix<N> J(model.size(),model.size()); // Jacobian matrix
      Vector<N> d_newton(model.size());              // solution of linear system, which returns the newton direction
      Vector<N> s(model.size());              // scaling factors
      Vector<size_type> p(model.size());                 // row permutations
      Vector<size_type> q(model.size());                 // column permutations

      model.F(x,F);                                     // compute nonlinear residual
     
      Real R0(std::abs(norm(F)));                          // norm of initial residual
      Real R(R0);                                // current residual norm

      if (verbosity>=1)
      {
        std::cout << "Newton Dog Leg Cauchy " 
                  << "   norm=" << std::scientific << std::showpoint 
                  << std::setprecision(4) << R0
                  << std::endl;
      }

      converged = false;

      //Initial trust radius
      trust_region = initial_trust_radius;

      //save the type of the direction which was used in each step
      std::string direction_type;

      for (size_type i=1; i<=maxit; i++)                // do iterations
        {
          iteration_points.push_back(x);
          norms.push_back(R);
          loss.push_back(0.5 * R * R );
          trust_radius_list.push_back(trust_region);

          //check convergence
          if (R<=abslimit)
            {
              if (verbosity>=1)
              {
                std::cout << "Newton Dog Leg Cuchy converged in "  << i << " steps"
                            << " reduction=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << R/R0
                            << std::endl;
              }
              iterations_taken = i;
              converged = true;
              if(filename.size()!=0){
                saveDatainFile(filename,iteration_points, newton_directions, steepest_descent_directions, dog_leg_directions, trust_radius_list, norms, reductions, loss);
              }
              return;
            } 

          model.F_x(x,J);                               // compute Jacobian matrix
          model.F(x,F);                                 // compute residual
          
          hdnum::Vector<N> jF = J.transpose() * F;      // compute J^T * F
          hdnum::DenseMatrix<N> B = J.transpose() * J;  // compute J^T * J

          hdnum::DenseMatrix<N> A = J;

          row_equilibrate(A,s);                         // equilibrate rows
          lr_fullpivot(A,p,q);                          // LR decomposition of A
          d_newton = N(0.0);                                   // clear solution
          apply_equilibrate(s,F);                       // equilibration of right hand side
          permute_forward(p,F);                         // permutation of right hand side
          solveL(A,F,F);                                // forward substitution
          solveR(A,d_newton,F);                                // backward substitution
          permute_backward(q,d_newton);                        // backward permutation

          d_newton *= -1.0;

          //step direction 
          Vector<N> d(model.size());

          //check if newton direction can be applied
          if(std::abs(norm(d_newton)) <= trust_region){
            d = d_newton;
            direction_type = "Newton direction";

            newton_directions.push_back(d_newton);
            Vector<N> zeros(model.size());
            steepest_descent_directions.push_back(zeros);
            dog_leg_directions.push_back(zeros);
          }
          else{
            Vector<N> zeros(model.size());
            newton_directions.push_back(zeros);
            double norm_jF = std::abs(norm(jF));
            Vector<N> product_J_jF= J.transpose() * jF;
            double norm_B_jF = std::abs(norm(product_J_jF));

            double alpha = norm_jF * norm_jF / (norm_B_jF * norm_B_jF);

            Vector<N> d_cauchy = jF;
            d_cauchy *= -alpha;
            
            //check if steepest descent direction can be applied
            if(std::abs(norm(d_cauchy)) >= trust_region){
              d_cauchy /= std::abs(norm(d_cauchy));
              d_cauchy *= trust_region;

              d = d_cauchy;
              direction_type = "Steepest descent direction";

              steepest_descent_directions.push_back(d_cauchy);
              Vector<N> zeros(model.size());
              dog_leg_directions.push_back(zeros);
            }
            else{ 
              Vector<N> zeros(model.size());
              steepest_descent_directions.push_back(zeros);
              //apply dog leg direction by soolving the quadratic equation ||d_cauchy + tau * (d_newton-d_cauchy)|| = trust_region
              Vector<N> dn_minus_dc = d_newton - d_cauchy;
  
              double a = std::abs(norm(dn_minus_dc)) * std::abs(norm(dn_minus_dc));
              N product_dc_dn_minus_dc = d_cauchy * dn_minus_dc;
              N b = 2.0 * product_dc_dn_minus_dc;

              N dc_dc = d_cauchy * d_cauchy;
              N c = -trust_region * trust_region + dc_dc;

              N root = sqrt(b*b - 4 * a * c);
              N solution = (-b + root) / (2 * a);

              Vector<N> product_dn_minus_dc_solution = dn_minus_dc;
              product_dn_minus_dc_solution *= solution;

              d =  d_cauchy + product_dn_minus_dc_solution;
              direction_type = "Dog leg direction";

              dog_leg_directions.push_back(d);
            } 
          }

          /*compute ratio of actual reduction and predicted reduction(How good is the approximated model?)
          * (f(x) - f(x+d)) / (m(0) - m(d))
          * with m(p) = f(x) + (J^T * F(x))^T * p + p^T * J^T*J * p
          */
          Vector<N> new_x = x + d;

          Vector<N> res_old(model.size());
          model.F(x,res_old);
          Vector<N> res_new(model.size());
          model.F(new_x, res_new);
          
          Real newR(std::abs(norm(res_new)));

          //reduction for the output
          Real red(1.0);

          N product_res_old_res_old = res_old * res_old;
          N product_res_new_res_new = res_new * res_new;
          N actual_reduction = 0.5 * product_res_old_res_old - 0.5 *  product_res_new_res_new;

          N product_jF_d = jF * d;
          Vector<N> product_B_d = B * d;
          N product_d_B_d = d * product_B_d;
          product_d_B_d *= 0.5;
          N predicted_reduction = - (product_jF_d + product_d_B_d);

          N ro = actual_reduction / predicted_reduction;

          if(std::abs(ro) < n2){
            trust_region = t1 * trust_region;
          }
          else if(std::abs(ro) > n3 && std::abs(norm(d)) == trust_region){
            trust_region = std::min(t2 * trust_region, max_radius);
          }
          if(std::abs(ro) > n1){
            x = new_x;
            red = newR/R;
            R = newR;
          }

          reductions.push_back(red);

          if (verbosity>=2)
          {
            std::cout << "  step"  << std::setw(3) << i 
                      << " norm=" << std::scientific << std::showpoint 
                      << std::setprecision(4) << R
                      << " red=" << std::scientific << std::showpoint 
                      << std::setprecision(4) << red
                      << "  " << direction_type
                      << std::endl;
          }
         

          if (R<=reduction*R0)
            {
              if (verbosity>=1)
              {
                std::cout << "Newton Dog Leg Cuchy converged in "  << i << " steps"
                            << " reduction=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << R/R0
                            << std::endl;
              }
              iterations_taken = i;
              converged = true;
              if(filename.size()!=0){
                saveDatainFile(filename,iteration_points, newton_directions, steepest_descent_directions, dog_leg_directions, trust_radius_list, norms, reductions, loss);
              }
              return;
            }

          if(i == maxit){
            iterations_taken = i;
            if (verbosity>=1){
              std::cout << "Newton Dogleg Cauchy not converged within " << maxit << " iterations" << std::endl;
              if(filename.size()!=0){
                saveDatainFile(filename,iteration_points, newton_directions, steepest_descent_directions, dog_leg_directions, trust_radius_list, norms, reductions, loss);
              }
            }
          }

        }

    }

        
    bool has_converged () const
    {
      return converged;
    }
    size_type iterations() const {
      return iterations_taken;
    }


    private:
      size_type maxit;

      mutable size_type iterations_taken = -1;
      size_type linesearchsteps;
      size_type verbosity;
      double reduction;
      double abslimit;
      mutable bool converged;

      double trust_region = 1.0;
      double initial_trust_radius = 1.0;
      double max_radius = 1.0;

      double n1 = 0.001;
      double n2 = 0.25;
      double n3 = 0.75;
      double t1 = 0.25;
      double t2 = 2.0;
  };

  class ProjectedNewton{
    typedef std::size_t size_type;

    template<typename N>
    void saveDatainFile(std::string filename, std::vector<std::vector<N>> iteration_points){
      std::fstream file(filename.c_str(),std::ios::out);
      for(int i = 0; i< iterations_taken;++i){
        Vector<N> x(2);
        x[0] = iteration_points[i][0]; 
        x[1] = iteration_points[i][1];

        file << x[0] << "   " << x[1] << "\n";
      }
    }

    public:
      
    ProjectedNewton ()
      : maxit(25), linesearchsteps(30), converged(false)
    {}

    void set_maxit (size_type n)
    {
      maxit = n;
    }

    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x, std::string filename="") 
    {
      typedef typename M::number_type N;

      std::vector<std::vector<N>> iteration_points;

      DenseMatrix<N> constraints = model.A;
      Vector<N> lowerbounds = model.lowerBounds;
      Vector<N> upperbounds = model.upperBounds;
      
      DenseMatrix<N> activeSet = constraints.sub(0,0, model.size(), model.size());
      Vector<N> activelowerbounds = lowerbounds.sub(0, model.size());
      Vector<N> activeupperbounds = upperbounds.sub(0, model.size());

      int u = 0;

      int k = 0;
      while(k <=5){
        Vector<N> y = activeSet * x;
        for(int k = 0; k<y.size(); ++k){
          if(y[k] < activelowerbounds[k]){
            y[k] = activelowerbounds[k];
          }
          else if(y[k] > activeupperbounds[k]){
            y[k] = activeupperbounds[k];
          }
        }
        DenseMatrix<N> A = activeSet;

        //compute A^T
        DenseMatrix<N> A_T = activeSet;
        A_T = A_T.transpose(); 
        
        //LR decomposition of A 
        Vector<size_t> p_(2);
        Vector<size_t> q_(2);         
        Vector<N> s_(model.size()); 
        row_equilibrate(A,s_);                         // equilibrate rows
        lr_fullpivot(A,p_,q_);                          // LR decomposition of A

        //LR decompisition of A^T
        Vector<size_t> p(2);
        Vector<size_t> q(2);
        Vector<N> d(model.size());              
        Vector<N> s(model.size());            
        row_equilibrate(A_T,s);                         // equilibrate rows
        lr_fullpivot(A_T,p,q);                          // LR decomposition of A
        
        Vector<N> y_temp = y;
        x = N(0.0);
        apply_equilibrate(s_,y_temp);                       // equilibration of right hand side
        permute_forward(p_,y_temp);                         // permutation of right hand side
        solveL(A,y_temp,y_temp);                                // forward substitution
        solveR(A,x,y_temp);                                // backward substitution
        permute_backward(q_,x);                        // backward permutation

        //run newton iterations on transformed space y = Ax
        for(int i = 0; i< 1000; ++i){ 
          ++u;
          iteration_points.push_back(x);
          DenseMatrix<N> B(2,2);
          model.H(x,B);

          Vector<N> loss(1);
          model.f(x,loss);

          //solve A^T * d = g(x)
          Vector<N> gradient(model.size());         
          model.g(x, gradient);
          d = N(0.0);                                   // clear solution
          apply_equilibrate(s,gradient);                       // equilibration of right hand side
          permute_forward(p,gradient);                         // permutation of right hand side
          solveL(A_T,gradient,gradient);                                // forward substitution
          solveR(A_T,d,gradient);                                // backward substitution
          permute_backward(q,d);                        // backward permutation

          //determin step size
          double alpha = 1.0;
          int k = 0;
          double difference = 0.0;

          //line search
          while(true){
            if(k==500){
              break;
            }
            Vector<N> j = y;
            j.update(-alpha, B * d);

            //model.project_to_feasible_space(j);
            for(int k = 0; k<j.size(); ++k){
              if(j[k] < activelowerbounds[k]){
                j[k] = activelowerbounds[k];
              }
              else if(j[k] > activeupperbounds[k]){
                j[k] = activeupperbounds[k];
              }
            }

            Vector<N> j_temp = j;
            Vector<N> x_new(2);
            apply_equilibrate(s_,j_temp);                       // equilibration of right hand side
            permute_forward(p_,j_temp);                         // permutation of right hand side
            solveL(A,j_temp,j_temp);                                // forward substitution
            solveR(A,x_new,j_temp);                                // backward substitution
            permute_backward(q_,x_new);                        // backward permutation

            Vector<N> new_loss(1);
            model.f(x_new, new_loss);

            if(new_loss[0] < loss[0]){
              difference = norm(y-j);
              y = j;
              x = x_new;
              break;
            }

            alpha = alpha * 0.2;

            k = k+1;
          }

          //check convergence
          if(difference < 1e-12){
            break;
          }
        }

      int violated = model.constraintViolated(x);
      if(violated == -1){
        std::cout<<"converged"<<std::endl;
        std::cout<<"result: " << x <<std::endl;
        converged = true;
        iterations_taken = u;
        if(filename.size()!=0)
          saveDatainFile(filename, iteration_points);
        break;
      }
      else{
        DenseMatrix<N> temp = activeSet.sub(1,0, model.size()-1, model.size());
        Vector<N> a(model.size());
        for(int p = 0; p<= model.size();++p){
          a[p] = constraints(violated, p);
        }
        temp.addNewRow(a);

        activeSet = temp;
    
        Vector<N> temp1 = activelowerbounds.sub(1, model.size()-1);
        temp1.push_back(lowerbounds[violated]);

        activelowerbounds = temp1;

        Vector<N> temp2 = activeupperbounds.sub(1, model.size()-1);
        temp2.push_back(upperbounds[violated]);

        activeupperbounds = temp2;
      }
      k=k+1;

      }
    }


        
    bool has_converged () const
    {
      return converged;
    }
    size_type iterations() const {
      return iterations_taken;
    }


    private:
      size_type maxit;

      mutable size_type iterations_taken = -1;
      size_type linesearchsteps;
      mutable bool converged;
  };

} // namespace hdnum


#endif
