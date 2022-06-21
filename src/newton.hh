// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_NEWTON_HH
#define HDNUM_NEWTON_HH
#include <random>
#include <algorithm>

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
    {
    }

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
    template<typename Matrix>
    void F_x (const Vec& x, Matrix& result) const
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
          for (int i=0; i<result.rowsize(); i++){
            auto g = (Fz[i]-Fx[i])/dz;
            //std::cout << typeid(g).name() << std::endl;
            result[i][j] = (Fz[i]-Fx[i])/dz;
          }
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
   * @brief A generic nonlinear minimization problem class that can be set up with a lambda(Objective) definining min F(x),
   *        a lambda(gradient) defining the gradient of F and constraints defined by lower/upper bounds and a transformation matrix A defining 
   *        lower bound <= Ax <= upper bound.
   * 
   * @tparam Objective F(x) 
   * @tparam Gradient gradient of F
   * @tparam Hessian hessian of F
   * @tparam Vec the type of the Vector
   */
  template<typename Objective, typename Gradient, typename Hessian, typename Vec>
  class GenericNonlinearMinimizationProblem_Constrained
  {
  public:
    /** \brief export size_type */
    typedef std::size_t size_type;

    /** \brief export number_type */
    typedef typename Vec::value_type number_type;

    GenericNonlinearMinimizationProblem_Constrained (const Objective& o, const Gradient& g, const Hessian& h, const DenseMatrix<double> constraints, const Vector<double> lb, const Vector<double> up, const  Vec& x_, number_type eps_ = 1e-7)
      : objective(o), gradient(g), hessian(h), s(x_.size()), eps(eps_), upperBounds(up), lowerBounds(lb), A(constraints)
    {
        // check if own hessian was provided or not(If not then the hessian weill be a zero matrix and you can check the first row)
        for(int i = 0; i< x_.size() ; ++i){
          if(hessian(x_)[0][i] != 0.0){
            usingOwnHessian = true;
            break;
          }
        }
    }

    Vector<number_type>& getLowerBounds(){
      return lowerBounds;
    }

    Vector<number_type>& getUpperBounds(){
      return upperBounds;
    }

    DenseMatrix<number_type>& getTransformationMatrix(){
      return A;
    }

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
      if(usingOwnHessian){
        result = hessian(x);
      }
      else{
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
    }
    
    /* the active set A* is a sub set of the transformation matrix A and will be determined by the following rules: 
      (1) All active constraints (lower bound = Ax or Ax = upper bound) are contained in A*
      (2) Rang(A*) = size of x
    */
    void determineActiveSet(Vec& x, DenseMatrix<number_type>& activeSet, Vector<number_type>& activelowerbounds, Vector<number_type>& activeupperbounds) const{
      auto y = A * x;
      // project to feasible space
      for(int k = 0; k<y.size(); ++k){
        if(y[k] < lowerBounds[k]){
          y[k] = lowerBounds[k];
        }
        else if(y[k] > upperBounds[k]){
          y[k] = upperBounds[k];
        }
      }

      Vector<int> activeIndices; // save indices of active constraints
      Vector<int> inActiveIndices; // " "                  constraints

      for(int i = 0; i<y.size(); ++i){
        if(y[i] == lowerBounds[i] || y[i] == upperBounds[i]) 
          activeIndices.push_back(i);
        else
          inActiveIndices.push_back(i);
      }
      
      // check if Rang(A*) = s otherwise and add some random elements of the inactive constraints to the active set
      if(activeIndices.size() < s){ 
        std::random_device rd;
        std::mt19937 randomGenerator(rd());

        std::shuffle(inActiveIndices.begin(), inActiveIndices.end(), randomGenerator); // shuffle the inactive set

        for(int i = 0; i<inActiveIndices.size(); ++i){
          if(activeIndices.size() == s){ // stop if Rang(A*) = s
            break;
          }
          else{
            activeIndices.push_back(inActiveIndices[i]);
          }
        }
      }

      // build the active set by using the active indices
      for(int i = 0; i< x.size();++i){
        // add trivial constraints if necessary
        if(i >= activeIndices.size()){ 
          activelowerbounds[i] = std::numeric_limits<int>::min();
          activeupperbounds[i] = std::numeric_limits<int>::max();
          activeSet[i][i - activeIndices.size()] = 1.0;
        }
        else{
          int index = activeIndices[i];
          activelowerbounds[i] = lowerBounds[index];
          activeupperbounds[i] = upperBounds[index];

          for(int j = 0; j < x.size(); ++j){
            activeSet[i][j] = A[index][j];
          }
        }
      }
    }
    
    private:
    Objective objective; 
    Gradient gradient;
    Hessian hessian;

    size_t s;
    typename Vec::value_type eps;

    DenseMatrix<number_type> A;
    Vector<number_type> upperBounds;
    Vector<number_type> lowerBounds;

    // checks if hessian is provided
    bool usingOwnHessian= false;
  
  };

  /**
   * @brief A function returning a minimization problem class with constraints
   * 
   * @tparam Objective F(x) 
   * @tparam Gradient gradient of F
   * @tparam Hessian hessian of F
   * @tparam X the type of the Vector
   * @param constraints transformation matrix A
   * @param lb lower bound
   * @param up upper bound
   * @return GenericNonlinearMinimizationProblem_Constrained<Objective,Gradient,X> 
   */
  template<typename Objective,typename Gradient,typename Hessian, typename X>
  GenericNonlinearMinimizationProblem_Constrained<Objective,Gradient, Hessian, X> getNonlinearMinimizationProblem_Constrained(const Objective& o, const Gradient& g, const Hessian& h, const DenseMatrix<double>& constraints, const Vector<double>& lb, const Vector<double>& up, const  X& x, typename X::value_type eps = 1e-7)
  {
    return GenericNonlinearMinimizationProblem_Constrained<Objective, Gradient, Hessian, X>(o,g,h,constraints, lb, up, x,eps);
  }

  /** @brief Solve nonlinear problem using a damped Newton method

      The Newton solver is parametrized by a model. The model also
      exports all relevant types for types.
 
  */
  class Newton
  {
  //protected:
  protected:
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


    /**
     * @brief Solver method. 
     * 
     * @tparam M 
     * @tparam Vec 
     * @param model problem class
     * @param x initial solution. For real/complex values x is of type Vector/CVector
     * @param filename save results in a file when needed
     */
    template<class M, typename Vec>
    void solve (const M& model, Vec& x, std::string filename = "") const
    {
      // check if we have a complex problem
      const auto xTemp = x;
      ComplexChecker complexChecker;
      constexpr bool isComplex = complexChecker.isComplexVector(xTemp);  
      
      // in complex case, we still need to use real valued numbers for residual norms etc.
      using Real = typename M::number_type;

      // type of vector and matrix depends on whether we have a complex problem or not
      typedef typename std::conditional<isComplex, CVector<Real>, Vector<Real>>::type vectortype;
      typedef typename std::conditional<isComplex, CDenseMatrix<Real>, DenseMatrix<Real>>::type matrixtype;

      // for saving the data
      Vector<vectortype> iteration_points;
      Vector<vectortype> directions;
      Vector<Real> stepSizes;
      Vector<Real> norms;
      Vector<Real> reductions;
      Vector<Real> losses;

      vectortype r(model.size());              // residual
      matrixtype A(model.size(),model.size()); // Jacobian matrix
      vectortype y(model.size());              // temporary solution in line search
      vectortype z(model.size());              // solution of linear system
      vectortype s(model.size());              // scaling factors
      Vector<size_type> p(model.size());                 // row permutations
      Vector<size_type> q(model.size());                 // column permutations

      model.F(x,r);                                     // compute nonlinear residualz

      Real R0(norm(r));                          // norm of initial residual
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
          iteration_points.push_back(x);
          norms.push_back(R);
          losses.push_back(0.5 * R * R);

          // check absolute size of residual
          if (R<=abslimit)
            {
              converged = true;
              if(filename.size()>0)
                saveDataInFile(filename,iteration_points,directions, stepSizes, norms, reductions, losses);
              return;
            } 

          // solve Jacobian system for update
          model.F_x(x,A);                               // compute Jacobian matrix

          matrixtype temp = A;
          row_equilibrate(A,s);                         // equilibrate rows
          lr_fullpivot(A,p,q);                          // LR decomposition of A
          z = vectortype(model.size());                                      // clear solution
          apply_equilibrate(s,r);                       // equilibration of right hand side
          permute_forward(p,r);                         // permutation of right hand side
          solveL(A,r,r);                                // forward substitution
          solveR(A,z,r);                                // backward substitution
          permute_backward(q,z);                        // backward permutation
          z *= -1.0;
          directions.push_back(z);

          // line search
          Real lambda(1.0);                      // start with lambda=1
          for (size_type k=0; k<linesearchsteps; k++)
            {
              y = x;                              
              y.update(lambda,z);                       // y = x+lambda*z
              model.F(y,r);                             // r = F(y)
              Real newR(norm(r));                // compute norm
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
                  stepSizes.push_back(lambda);
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
                    if(filename.size()>0)
                      saveDataInFile(filename,iteration_points,directions, stepSizes, norms, reductions, losses);
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
              if(filename.size()>0)
                saveDataInFile(filename,iteration_points,directions, stepSizes, norms, reductions, losses);
              return;
            }
          if (i==maxit)
            {
              iterations_taken = i;
              if (verbosity>=1){
                std::cout << "Newton not converged within " << maxit << " iterations" << std::endl;
                if(filename.size()>0)
                  saveDataInFile(filename,iteration_points,directions, stepSizes, norms, reductions, losses);
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

  protected:
    size_type maxit;
    mutable size_type iterations_taken = -1;
    size_type linesearchsteps;
    size_type verbosity;
    double reduction;
    double abslimit;
    mutable bool converged;

    // save results in a file 
    template<typename vectortype, typename Real>
    void saveDataInFile(const std::string filename,const Vector<vectortype>& iterationPoints,const Vector<vectortype>& directions, const Vector<Real>& stepSizes,const Vector<Real>& norms, const Vector<Real>& reductions, const Vector<Real>& losses) const{
        std::fstream file(filename.c_str(),std::ios::out);
        file<< "#This file contains the outputs from the newton-raphson method solver. The outputs are ordered in the following way: \n";
        file<< "#Iteration point | direction | step size | norm | reduction | loss\n";

        // go through iterations
        for(int i = 0; i< iterations_taken;++i){
            for(auto x: iterationPoints[i])
              file << x << "   ";
            for(auto d: directions[i])
              file << d << "   ";

            file << stepSizes[i] << "   " << norms[i] << "   " << reductions[i]<< "   " << losses[i] << "\n"; 
        }
    }
  };

/** @brief  Solve nonlinear problem F(x) = 0 by defining an eqivalent minimizing problem min_x f(x) = 0.5 * F(x)^T * F(x) 
 * (Uses Trust region method where a subproblem is solved by different direction choices:
 *  Newton direction
 *  Cauchy point
 *  Dog leg direction)
 */
  class NewtonDogLegCauchy: public Newton{
    public:
      
    NewtonDogLegCauchy ()
      :maxTrustRadius(1.0), initialTrustRadius(1.0)
    {
      maxit = 25;
      verbosity = 0;
      reduction = 1e-14;
      abslimit = 1e-30;
      converged = false;
    }

    void setMaxTrustRadius(double max){
      maxTrustRadius = max;
    }

    void setInitialTrustRadius(double initial){
      if(initial > maxTrustRadius)
        initialTrustRadius = maxTrustRadius;
      else
        initialTrustRadius = initial;
    }

    /**
     * @brief Solver method. 
     * 
     * @tparam M 
     * @tparam Vec 
     * @param model problem class
     * @param x initial solution. For real/complex values x is of type Vector/CVector
     * @param filename save results in a file when needed
     */
    template<class M, typename Vec>
    void solve (const M& model, Vec& x, std::string filename="") const
    {
      // check if we have a complex problem
      const auto xTemp = x;
      ComplexChecker complexChecker;
      constexpr bool isComplex = complexChecker.isComplexVector(xTemp);  

      // in complex case, we still need to use real valued numbers for residual norms etc.
      using Real = typename M::number_type;

      // type of vector and matrix depends on whether we have a complex problem or not
      typedef typename std::conditional<isComplex, CVector<Real>, Vector<Real>>::type vectortype;
      typedef typename std::conditional<isComplex, CDenseMatrix<Real>, DenseMatrix<Real>>::type matrixtype;
      // for complex problems we also need complex values
      typedef typename std::conditional<isComplex, std::complex<Real>, Real>::type N;

      // for saving the results
      Vector<vectortype> iterationPoints;
      Vector<vectortype> newtonDirections;
      Vector<vectortype> steepestDescentDirections;
      Vector<vectortype> doglegDirections;
      Vector<double> trustRadiusList;
      Vector<Real> norms;
      Vector<Real> reductions;
      Vector<Real> losses;

      vectortype F(model.size());              // residual
      matrixtype J(model.size(),model.size()); // Jacobian matrix
      vectortype d_newton(model.size());              // solution of linear system, which returns the newton direction
      vectortype s(model.size());              // scaling factors
      Vector<size_type> p(model.size());                 // row permutations
      Vector<size_type> q(model.size());                 // column permutations

      model.F(x,F);                                     // compute nonlinear residual
     
      Real R0(norm(F));                          // norm of initial residual
      Real R(R0);                                // current residual norm
      
      if (verbosity>=1)
      {
        std::cout << "Newton Dog Leg Cauchy " 
                  << "   norm=" << std::scientific << std::showpoint 
                  << std::setprecision(4) << R0
                  << std::endl;
      }

      converged = false;

      double trustRadius = initialTrustRadius; // initial trust radius

      std::string direction_type; // save the type of the direction which was used in each iteration
      
      for (size_type i=1; i<=maxit; i++)                // do iterations
        {
          iterationPoints.push_back(x);
          norms.push_back(R);
          losses.push_back(0.5 * R * R );
          trustRadiusList.push_back(trustRadius);

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
              if(filename.size()!=0)
                saveDatainFile(filename,iterationPoints, newtonDirections, steepestDescentDirections, doglegDirections, trustRadiusList, norms, reductions, losses);
              return;
            } 
          
          /* Each iterations solves the trust region subproblem 
          *  min m(d) = f(x) + J^T * F * d + 0.5 * p^T * J^T * J * p s.t ||p|| <= trustRadius 
          *  by using different directions:
          *  - Newton direction: solves J^T d = - F(x) 
          *    => will be applied if direction is inside trust region
          *  - Cauchy point: takes the negative gradient with scaling factor alpha to minimize m(alpha *d) 
          *    => will be applied if direction is outside trust region
          *  - otherwise dogleg direction: interpolation of the endpoints from Newton and Cauchy direction
          */
          model.F_x(x,J);                               // compute Jacobian matrix
          model.F(x,F);                                 // compute residual

          matrixtype JTranspose =  getTranspose(J);
          vectortype jF = JTranspose * F;      // compute J^T * F (gradient of f)
          matrixtype B = JTranspose * J;  // compute J^T * J (Approximation of the Hessian of f)

          // copy J, otherwise it will be modified by the LR decomposition
          matrixtype A = J;

          // solve J^T d = - F(x) to get newton direction
          row_equilibrate(A,s);                         // equilibrate rows
          lr_fullpivot(A,p,q);                          // LR decomposition of A
          d_newton = vectortype(model.size());                                   // clear solution
          apply_equilibrate(s,F);                       // equilibration of right hand side
          permute_forward(p,F);                         // permutation of right hand side
          solveL(A,F,F);                                // forward substitution
          solveR(A,d_newton,F);                                // backward substitution
          permute_backward(q,d_newton);                        // backward permutation

          d_newton *= -1.0;

          vectortype d(model.size()); // direction to apply
          vectortype zeros(model.size()); // for directions which are not applied in current iteration

          // check if newton direction is inside trust region
          if(norm(d_newton) <= trustRadius){
            d = d_newton;

            direction_type = "Newton direction";
            newtonDirections.push_back(d_newton);
            steepestDescentDirections.push_back(zeros);
            doglegDirections.push_back(zeros);
          }
          else{
            newtonDirections.push_back(zeros);

            /* the cauchy point is given as d = - alpha * J^T F (steepest descent direction) and 
            * alpha = min( trustRadius / abs(norm(J^T F))) , abs(norm(J^T F)))**2 / ( (J^T F)^T * B * J^T F) )
            * For the implementation, instead of directly taking the min, first compute the second term and check if alpha * d
            * is outside the trust region or on the edge.
            */
            Real norm_jF = norm(jF);

            // (J^T*F)^T * B * J^T*F) = ||J^T*J^T*F||^2
            vectortype product_J_jF= JTranspose * jF;
            Real norm_B_jF = norm(product_J_jF);
            
            Real alpha_1 = norm_jF * norm_jF / (norm_B_jF * norm_B_jF);
            Real alpha_2 = trustRadius / norm_jF;
            Real alpha = std::min(alpha_1, alpha_2);

            vectortype d_cauchy = jF;
            d_cauchy *= -alpha;

            //check if cauchy point is outside trust region or on the edge
            if(norm(d_cauchy) >= trustRadius){
              d = d_cauchy;
              direction_type = "Steepest descent direction";

              steepestDescentDirections.push_back(d_cauchy);
              doglegDirections.push_back(zeros);
            }
            else{ 
              steepestDescentDirections.push_back(zeros);
              // apply dog leg direction by solving the quadratic equation ||d_cauchy + tau * (d_newton-d_cauchy)|| = trust_region
              vectortype dn_minus_dc = d_newton - d_cauchy;
  
              Real a = norm(dn_minus_dc);
              a*= a;
              N  product_dc_dn_minus_dc = d_cauchy * dn_minus_dc;
              N b = 2.0 * product_dc_dn_minus_dc;

              Real dc_dc = norm(d_cauchy);
              dc_dc *= dc_dc;
              Real c = -trustRadius * trustRadius + dc_dc;

              N root = sqrt(b*b - 4 * a * c);
              N solution = (-b + root) / (2 * a);

              vectortype product_dn_minus_dc_solution = dn_minus_dc;
              product_dn_minus_dc_solution *= solution;

              d =  d_cauchy + product_dn_minus_dc_solution;
              direction_type = "Dog leg direction";

              doglegDirections.push_back(d);
            } 
          }

          /*compute ratio of actual reduction and predicted reduction(How good is the approximated model?):
          * actualReduction/predicted Reduction = (f(x) - f(x+d)) / (m(0) - m(d)), 
          * with m(d) = f(x) + (J^T * F(x))^T * d + d^T * B * d,  B = J^T * J and f(x) = 0.5 * F(x)^T * F(x) 
          */
          vectortype new_x = x + d; // perform update
          
          vectortype res_old(model.size());
          model.F(x,res_old); // F(x)
          vectortype res_new(model.size());
          model.F(new_x, res_new); // F(x+d)
          
          Real newR(norm(res_new)); // norm of new residual
          Real red(1.0); // reduction for the output

          Real actual_reduction = 0.5 * R * R - 0.5 *  newR * newR; //  (f(x) - f(x+d))

          Real product_jF_d = std::real(jF * d); // (J^T * F(x))^T * d 
          vectortype product_B_d = B * d;
          
          Real product_d_B_d= std::real(d * product_B_d);
          product_d_B_d *= 0.5; // 0.5 * d^T * B * d

          Real predicted_reduction = - (product_jF_d + product_d_B_d); 
          Real ro = actual_reduction / predicted_reduction;
          Real absro = std::abs(ro);

          if(absro < n2){
            trustRadius = t1 * trustRadius; // bad approximation
          }
          else if(absro > n3 && norm(d) == trustRadius){
            trustRadius = std::min(t2 * trustRadius, maxTrustRadius); // good approximation
          }

          // check if approximation is good enough to perform update
          if(absro > n1){
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
                std::cout << "Newton Dog Leg Cauchy converged in "  << i << " steps"
                            << " reduction=" << std::scientific << std::showpoint 
                            << std::setprecision(4) << R/R0
                            << std::endl;
              }
              iterations_taken = i;
              converged = true;
              if(filename.size()!=0)
                saveDatainFile(filename,iterationPoints, newtonDirections, steepestDescentDirections, doglegDirections, trustRadiusList, norms, reductions, losses);
              return;
            }

          if(i == maxit){
            iterations_taken = i;
            if (verbosity>=1){
              std::cout << "Newton Dogleg Cauchy not converged within " << maxit << " iterations" << std::endl;
              if(filename.size()!=0){}
                saveDatainFile(filename,iterationPoints, newtonDirections, steepestDescentDirections, doglegDirections, trustRadiusList, norms, reductions, losses);
            }
          }

        }
    }

    private:
      double initialTrustRadius;
      double maxTrustRadius;

      // approximation values
      double n1 = 0.001;
      double n2 = 0.25;
      double n3 = 0.75;

      // scale factors for trust radius
      double t1 = 0.25;
      double t2 = 2.0;

      // save results in a file 
      template<typename vectortype, typename Real>
      void saveDatainFile(const std::string& filename, const Vector<vectortype>& iterationPoints, const Vector<vectortype>& newtonDirections, const Vector<vectortype>& steepestDescentDirections,const Vector<vectortype>& doglegDirections,const Vector<double>& trustRadiusList,const Vector<Real>& norms,const Vector<Real>& reductions,const Vector<Real>& losses) const{
        std::fstream file(filename.c_str(),std::ios::out);
        file<< "#This file contains the outputs from the newton-dog-leg-cauchy method solver. The outputs are ordered in the following way: \n";
        file<< "#Iteration point | directions(Newton-Steepest-Dog leg) | trust radius | norm | reduction | loss \n";
        
        // go through iterations
        for(int i = 0; i< iterations_taken;++i){
          for(auto x: iterationPoints[i])
            file << x << "   ";
          for(auto dNewton: newtonDirections[i])
            file << dNewton << "   ";
          for(auto dSteepest: steepestDescentDirections[i])
            file << dSteepest << "   ";
          for(auto dDogleg: doglegDirections[i])
            file << dDogleg << "   ";

          file << trustRadiusList[i] << "   " << norms[i] << "   " << reductions[i] << "   " << losses[i] <<"\n";
        }
      }

      template<typename N>
      DenseMatrix<N> getTranspose(const DenseMatrix<N>& A) const{
        return A.transpose();
      }

      //For complex case we use the conjugate transposed matrix
      template<typename N>
      CDenseMatrix<N> getTranspose(const CDenseMatrix<N>& A) const{
        return A.ctranspose();
      }

  };


/**
 * @brief Solve nonliner minimization problem min f(x) s.t. b1 <= Ax <= b2 with the projected newton method by performing the updates
 *        on the transformed variable y = Ax. A must have Rang(A) = size of x.
 * 
 */
  class ProjectedNewton : public Newton{
    typedef std::size_t size_type;

    public:
    ProjectedNewton ()
    {
      maxit = 250;
      linesearchsteps = 500;
      converged = false;
    }

    
    // Solver method. The last parameter is optional to save the results in a file
    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x, std::string filename="") const
    {
      typedef typename M::number_type N;
      using Real = typename std::conditional<std::is_same<std::complex<double>, N>::value, double, N>::type;

      Vector<Vector<N>> iteration_points; // for saving the results

      DenseMatrix<N> activeSet(model.size(), model.size()); // active set A*
      Vector<N> activelowerbounds(model.size());
      Vector<N> activeupperbounds(model.size());

      for(int i=0; i<=maxit; ++i){
        iteration_points.push_back(x);
        model.determineActiveSet(x, activeSet, activelowerbounds, activeupperbounds);

        // compute y = A* x
        Vector<N> y = activeSet * x;
        // project to feasible space
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
        
        // LR decomposition of A 
        Vector<size_t> p_(model.size());
        Vector<size_t> q_(model.size());         
        Vector<N> s_(model.size()); 
        row_equilibrate(A,s_);                         // equilibrate rows
        lr_fullpivot(A,p_,q_);                          // LR decomposition of A

        // LR decompisition of A^T
        Vector<size_t> p(model.size());
        Vector<size_t> q(model.size());
        Vector<N> d(model.size());              
        Vector<N> s(model.size());            
        row_equilibrate(A_T,s);                         // equilibrate rows
        lr_fullpivot(A_T,p,q);                          // LR decomposition of A
        
        // get x again to compute loss and Hessian
        Vector<N> y_temp = y;
        x = N(0.0);
        apply_equilibrate(s_,y_temp);                       // equilibration of right hand side
        permute_forward(p_,y_temp);                         // permutation of right hand side
        solveL(A,y_temp,y_temp);                                // forward substitution
        solveR(A,x,y_temp);                                // backward substitution
        permute_backward(q_,x);                        // backward permutation

        DenseMatrix<N> H(model.size(),model.size()); // Hessian
        model.H(x,H);

        Vector<N> loss(1); // loss(x) = f(x)
        model.f(x,loss);

        // solve z = H^-1 * (A*^T)^-1 * g(x) by (1) A*^T * d = g(x) and (2) H * z = d 

        //solve A*^T * d = g(x)
        Vector<N> gradient(model.size());         
        model.g(x, gradient);

        d = N(0.0);                                   // clear solution
        apply_equilibrate(s,gradient);                       // equilibration of right hand side
        permute_forward(p,gradient);                         // permutation of right hand side
        solveL(A_T,gradient,gradient);                                // forward substitution
        solveR(A_T,d,gradient);                                // backward substitution
        permute_backward(q,d);                        // backward permutation
        
        //solve H * z = d 
        Vector<N> z(model.size());

        // LR decomposition of H 
        Vector<size_t> p_H(model.size());
        Vector<size_t> q_H(model.size());         
        Vector<N> s_H(model.size()); 
        row_equilibrate(H,s_H);                         // equilibrate rows
        lr_fullpivot(H,p_H,q_H);                          // LR decomposition of A
        apply_equilibrate(s_H,d);                       // equilibration of right hand side
        permute_forward(p_H,d);                         // permutation of right hand side
        solveL(H,d,d);                                // forward substitution
        solveR(H,z,d);                                // backward substitution
        permute_backward(q_H,z);                        // backward permutation

        double alpha = 1e10; // step size
        int k = 0;
        Real difference = 0.0; // save difference between y and new updated y
        bool reduced = false; // check if line search was successfull

        while(true){
          if(k==linesearchsteps){
            break;
          }
          Vector<N> yNew = y;
          yNew.update(-alpha, z); // perform update y = y - alpha + z,  z = H^-1 * (A*^T)^-1 * g(x)

          // project to feasible space
          for(int k = 0; k<yNew.size(); ++k){
            if(yNew[k] < activelowerbounds[k]){
              yNew[k] = activelowerbounds[k];
            }
            else if(yNew[k] > activeupperbounds[k]){
              yNew[k] = activeupperbounds[k];
            }
          }
          // compute x
          Vector<N> yNewTemp = yNew;
          Vector<N> xNew(model.size());
          apply_equilibrate(s_,yNewTemp);                       // equilibration of right hand side
          permute_forward(p_,yNewTemp);                         // permutation of right hand side
          solveL(A,yNewTemp,yNewTemp);                                // forward substitution
          solveR(A,xNew,yNewTemp);                                // backward substitution
          permute_backward(q_,xNew);                        // backward permutation

          Vector<N> newLoss(1);
          model.f(xNew, newLoss);

          // check if there was a reduction
          if(newLoss[0] < loss[0]){
            difference = norm(y-yNew);
            y = yNew;
            x = xNew;
            reduced = true;
            break;
          }

          alpha = alpha * 0.5; // decrease step size
          k = k+1;
        }
        // check for convergence
        if(difference < 1e-14 && reduced){
          converged = true;
          iterations_taken = i;
          if(filename.size()!=0)
            saveDatainFile(filename, iteration_points);
          return;
        }

        if(i == maxit){
          std::cout<<"ee"<<std::endl;
          iterations_taken = i;
          if(filename.size()!=0)
            saveDatainFile(filename, iteration_points);
          return;
        }
      }
    }

    private:
    // save results in a file 
    template<typename N>
    void saveDatainFile(const std::string& filename,const Vector<Vector<N>>& iteration_points) const {
      std::fstream file(filename.c_str(),std::ios::out);

      // go through iterations
      for(int i = 0; i< iterations_taken;++i){
        for(auto x: iteration_points[i])
          file << x << "   ";
        file << "\n";
      }
    }
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


} // namespace hdnum


#endif
