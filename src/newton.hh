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


  /** @brief Solve nonlinear problem using a damped Newton method

      The Newton solver is parametrized by a model. The model also
      exports all relevant types for types.

  */
  class Newton
  {
    typedef std::size_t size_type;

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
      
    void save_data_in_file(std::string filename){
        std::fstream file(filename.c_str(),std::ios::out);
        file<< "#This file contains the outputs from the newton-raphson method solver. The outputs are ordered in the following way: \n";
        file<< "#Iteration point (X,Y)) | directions | step size | norm | reduction \n";
        for(int i = 0; i< iterations_taken;++i){
            Vector<double> x(2);
            x[0] = iteration_points[i][0]; 
            x[1] = iteration_points[i][1];

            Vector<double> direction(2);
            direction[0] = directions[i][0];
            direction[1] = directions[i][1];
            file << x[0] << "   " << x[1] << "   " << direction[0]<< "   " << direction[1] << "   "<< step_sizes[i] << "   " << norms[i] << "   " << reductions[i]<< "   " << loss[i] << "\n"; 

        }
    }

    std::vector<std::vector<double>> get_iteration_points(){
      return iteration_points;
    }

    std::vector<std::vector<double>> get_directions(){
      return directions;
    }

    std::vector<double> get_norms(){
      return norms;
    }

    std::vector<double> get_reductions(){
      return reductions;
    }

    void clear_data(){
      iteration_points.clear();
      directions.clear();
      step_sizes.clear();
      norms.clear();
      reductions.clear();
      loss.clear();
    }
    //! do one step
    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x) 
    {
      clear_data();

      typedef typename M::number_type N;
      // In complex case, we still need to use real valued numbers for residual norms etc.
      using Real = typename std::conditional<std::is_same<std::complex<double>, N>::value, double, N>::type;
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
          norms.push_back(norm(r));
          loss.push_back(0.5 * (r*r));

          // check absolute size of residual
          if (R<=abslimit)
            {
              converged = true;
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
                  if (verbosity>=3)
                    std::cout << "    line search not converged within " << linesearchsteps << " steps" << std::endl;
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
              return;
            }
          if (i==maxit)
            {
              iterations_taken = i;
              if (verbosity>=1)
                std::cout << "Newton not converged within " << maxit << " iterations" << std::endl;
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

    std::vector<std::vector<double>> iteration_points;
    std::vector<std::vector<double>> directions;
    std::vector<double> step_sizes;
    std::vector<double> norms;
    std::vector<double> reductions;
    std::vector<double> loss;

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
          std::cout<<"i: " << i<< std::endl;
          std::cout<<"x: " << x << std::endl;
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

/** @brief  Solver nonlinear problem F(x) = 0 by defining an eqivalent minimizing problem arg min_x f(x) = 0.5 * F(x)^T * F(x).
 * (Uses Trust region method where a subproblem is solved by different direction choices:
 *  Newton direction
 *  Steepest descent direction
 *  Dog leg direction)
 */
  class Newton_Dog_Leg_Cauchy{
    typedef std::size_t size_type;
    public:
      
    Newton_Dog_Leg_Cauchy ()
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

    void save_data_in_file(std::string filename){
      std::fstream file(filename.c_str(),std::ios::out);
      file<< "#This file contains the outputs from the newton-dog-leg-cauchy method solver. The outputs are ordered in the following way: \n";
      file<< "#Iteration point (X,Y)) | directions(Newton-Steepest-Dog leg) | trust radius | norm | reduction \n";
      for(int i = 0; i< iterations_taken;++i){
        Vector<double> x(2);
        x[0] = iteration_points[i][0]; 
        x[1] = iteration_points[i][1];

        Vector<double> newton_direction(2);
        newton_direction[0] = newton_directions[i][0];
        newton_direction[1] = newton_directions[i][1];

        Vector<double> steepest_descent_direction(2);
        steepest_descent_direction[0] = steepest_descent_directions[i][0];
        steepest_descent_direction[1] = steepest_descent_directions[i][1];

        Vector<double> dog_leg_direction(2);
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

    void clear_data(){
      iteration_points.clear();
      newton_directions.clear();
      steepest_descent_directions.clear();
      dog_leg_directions.clear();

      trust_radius_list.clear();
      norms.clear();
      reductions.clear();
      loss.clear();
    }

    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x) 
    {
      clear_data();

      typedef typename M::number_type N;

      // In complex case, we still need to use real valued numbers for residual norms etc.
      using Real = typename std::conditional<std::is_same<std::complex<double>, N>::value, double, N>::type;
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
          loss.push_back(0.5 * (R*R));
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
          if(norm(d_newton) <= trust_region){
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
            N norm_jF = norm(jF);
            Vector<N> product_J_jF= J.transpose() * jF;
            N norm_B_jF = norm(product_J_jF);

            N alpha = norm_jF * norm_jF / (norm_B_jF * norm_B_jF);

            Vector<N> d_cauchy = jF;
            d_cauchy *= -alpha;
            
            //check if steepest descent direction can be applied
            if(norm(d_cauchy) >= trust_region){
              d_cauchy /= norm(d_cauchy);
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
  
              N a = dn_minus_dc * dn_minus_dc;
              N product_dc_dn_minus_dc = d_cauchy * dn_minus_dc;
              N b = 2 * product_dc_dn_minus_dc;

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
          N red = 1.0;

          N product_res_old_res_old = res_old * res_old;
          N product_res_new_res_new = res_new * res_new;
          N actual_reduction = 0.5 * product_res_old_res_old - 0.5 *  product_res_new_res_new;

          N product_jF_d = jF * d;
          Vector<N> product_B_d = B * d;
          N product_d_B_d = d * product_B_d;
          product_d_B_d *= 0.5;
          N predicted_reduction = - (product_jF_d + product_d_B_d);

          N ro = actual_reduction / predicted_reduction;

          if(ro < n2){
            trust_region = t1 * trust_region;
          }
          else if(ro > n3 && norm(d) == trust_region){
            trust_region = std::min(t2 * trust_region, max_radius);
          }
          if(ro > n1){
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
              return;
            }

          if(i == maxit){
            iterations_taken = i;
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


      std::vector<std::vector<double>> iteration_points;
      std::vector<std::vector<double>> newton_directions;
      std::vector<std::vector<double>> steepest_descent_directions;
      std::vector<std::vector<double>> dog_leg_directions;

      std::vector<double> trust_radius_list;
      std::vector<double> norms;
      std::vector<double> reductions;
      std::vector<double> loss;
  };

} // namespace hdnum





#endif
