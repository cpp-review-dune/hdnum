// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_RUNGEKUTTA_HH
#define HDNUM_RUNGEKUTTA_HH

#include "vector.hh"
#include "newton.hh"

/** @file
 *  @general Runge-Kutta solver
 */

namespace hdnum {
  /** @brief Nonlinear problem we need to solve to do one step of an implicit Runge Kutta method
   */
  template<class M>
  class ImplicitRungeKuttaStepProblem
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores parameter lambda
    ImplicitRungeKuttaStepProblem (const M& model_,
                                   DenseMatrix<number_type> A_,
                                   Vector<number_type> b_,
                                   Vector<number_type> c_,
                                   time_type t_,
                                   Vector<number_type> u_,
                                   time_type dt_)
        : model(model_) , u(model.size())
      {
        A = A_;
        b = b_;
        c = c_;
        s = A_.rowsize ();
        dt = dt_;
        n = model.size();
        t = t_;
        u = u_;
      }

    //! return number of componentes for the model
    std::size_t size () const
    {
      return n*s;
    }

    //! model evaluation
    void F (const Vector<number_type>& x, Vector<number_type>& result) const
    {
      Vector<Vector<number_type> > xx (s);
      for (int i = 0; i < s; i++)
      {
        xx[i].resize(n,number_type(0));
        for(int k = 0; k < n; k++)
        {
          xx[i][k] = x[i*n + k];
        }
      }
      Vector<Vector<number_type> > f (s);
      for (int i = 0; i < s; i++)
      {
        f[i].resize(n, number_type(0));
        model.f(t + c[i] * dt, u + xx[i], f[i]);
      }
      Vector<Vector<number_type> > hr (s);
      for (int i = 0; i < s; i++)
      {
        hr[i].resize(n, number_type(0));
      }
      for (int i = 0; i < s; i++)
      {
        Vector<number_type> sum (n, number_type(0));
        for (int j = 0; j < s; j++)
        {
          sum.update(dt*A[i][j], f[j]);
        }
        hr[i]  = xx[i] - sum;
      }
      //translating hr into result
      for (int i = 0; i < s; i++)
      {
        for (int j = 0; j < n; j++)
        {
          result[i*n + j] = hr[i][j];
        }
      }
    }

    //! jacobian evaluation needed for newton in implicite solvers
    void F_x (const Vector<number_type>& x, DenseMatrix<number_type>& result) const
    {
      Vector<Vector<number_type> > xx (s);
      for (int i = 0; i < s; i++)
      {
        xx[i].resize(n);
        for(int k = 0; k < n; k++)
        {
          xx[i][k] = x[i*n + k];
        }
      }
      DenseMatrix<number_type> I (n, n, 0.0);
      for (int i = 0; i < n; i++)
      {
        I[i][i] = 1.0;
      }
      for (int i = 0; i < s; i++)
      {
        for (int j = 0; j < s; j++)
        {
          DenseMatrix<number_type> J (n, n, number_type(0));
          DenseMatrix<number_type> H (n, n, number_type(0));
          model.f_x(t+c[j]*dt, u + xx[j],H);
          J.update(-dt*A[i][j],H);
          if(i==j)                                //add I on diagonal
          {
            J+=I;
          }
          for (int k = 0; k < n; k++)
          {
            for (int l = 0; l < n; l++)
            {
              result[n * i + k][n * j + l] = J[k][l];
            }
          }
        }
      }
    }

  private:
    const M& model;
    time_type t, dt;
    Vector<number_type> u;
    int n, s;						    // dimension of matrix A and model.size
    DenseMatrix<number_type> A;				// A, b, c as in the butcher tableau
    Vector<number_type> b;
    Vector<number_type> c;
  };


  /** @brief classical Runge-Kutta method (order n with n stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M The model type
      \tparam S (Nonlinear) solver (default is Newton)
  */
  template<class M, class S = Newton>
  class RungeKutta
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    RungeKutta (const M& model_,
                DenseMatrix<number_type> A_,
                Vector<number_type> b_,
                Vector<number_type> c_)
      : model(model_), u(model.size()), w(model.size()), K(A_.rowsize ())
    {
      A = A_;
      b = b_;
      c = c_;
      s = A_.rowsize ();
      n = model.size();
      model.initialize(t,u);
      dt = 0.1;
      for (int i = 0; i < s; i++)
      {
        K[i].resize(n, number_type(0));
      }
      sigma = 0.01;
      verbosity = 0;

      if (A_.rowsize()!=A_.colsize())
      HDNUM_ERROR("need square and nonempty matrix");
      if (A_.rowsize()!=b_.size())
      HDNUM_ERROR("vector incompatible with matrix");
      if (A_.colsize()!=c_.size())
      HDNUM_ERROR("vector incompatible with matrix");
   }

   //! constructor stores reference to the model
   RungeKutta (const M& model_,
               DenseMatrix<number_type> A_,
               Vector<number_type> b_,
               Vector<number_type> c_,
               number_type sigma_)
     : model(model_), u(model.size()), w(model.size()), K(A_.rowsize ())
   {
     A = A_;
     b = b_;
     c = c_;
     s = A_.rowsize ();
     n = model.size();
     model.initialize(t,u);
     dt = 0.1;
     for (int i = 0; i < s; i++)
     {
       K[i].resize(n, number_type(0));
     }
     sigma = sigma_;
     verbosity = 0;
     if (A_.rowsize()!=A_.colsize())
     HDNUM_ERROR("need square and nonempty matrix");
     if (A_.rowsize()!=b_.size())
     HDNUM_ERROR("vector incompatible with matrix");
     if (A_.colsize()!=c_.size())
     HDNUM_ERROR("vector incompatible with matrix");
  }

  //! set time step for subsequent steps
  void set_dt (time_type dt_)
  {
    dt = dt_;
  }

  //! test if method is explicit
  bool check_explicit ()
  {
    bool is_explicit = true;
    for (int i = 0; i < s; i++)
    {
      for (int j = i; j < s; j++)
      {
        if (A[i][j] != 0.0)
        {
          is_explicit = false;
        }
      }
    }
    return is_explicit;
  }

  //! do one step
  void step ()
  {
    if (check_explicit())
    {
      // compute k_1
      w = u;
      model.f(t, w, K[0]);
      for (int i = 0; i < s; i++)
      {
        Vector<number_type> sum (K[0].size(), 0.0);
        sum.update(b[0], K[0]);
        //compute k_i
        for (int j = 0; j < i+1; j++)
        {
          sum.update(A[i][j],K[j]);
        }
        Vector<number_type> wert = w.update(dt,sum);
        model.f(t + c[i]*dt, wert, K[i]);
        u.update(dt *b[i], K[i]);
      }
    }
    if (not check_explicit())
    {
      // In the implicit case we need to solve a nonlinear problem
      // to do a time step.
      ImplicitRungeKuttaStepProblem<M> problem(model, A, b, c, t, u, dt);
      bool last_row_eq_b = true;
      for (int i = 0; i<s; i++)
      {
        if (A[s-1][i] != b[i])
        {
          last_row_eq_b = false;
        }
      }

      // Solve nonlinear problem and determine coefficients
      S solver;
      solver.set_maxit(2000);
      solver.set_verbosity(verbosity);
      solver.set_reduction(1e-10);
      solver.set_abslimit(1e-10);
      solver.set_linesearchsteps(10);
      solver.set_sigma(0.01);
      Vector<number_type> zij (s*n,0.0);
      solver.solve(problem,zij);


      DenseMatrix<number_type> Ainv (s,s,number_type(0));
      if (not last_row_eq_b)
      {
        // Compute LR decomposition of A
        Vector<number_type> w (s, number_type(0));
        Vector<number_type> x (s, number_type(0));
        Vector<number_type> z (s, number_type(0));
        Vector<std::size_t> p(s);
        Vector<std::size_t> q(s);
        DenseMatrix<number_type> Temp (s,s,0.0);
        Temp = A;
        row_equilibrate(Temp,w);
        lr_fullpivot(Temp,p,q);

        // Use LR decomposition to calculate inverse of A
        for (int i=0; i<s; i++)
        {
          Vector<number_type> e (s, number_type(0));
          e[i]=number_type(1);
          apply_equilibrate(w,e);
          permute_forward(p,e);
          solveL(Temp,e,e);
          solveR(Temp,z,e);
          permute_backward(q,z);
          for (int j = 0; j < s; j++)
          {
	        Ainv[j][i] = z[j];
          }
        }
      }

      Vector<Vector<number_type> > Z (s, 0.0);
      for(int i=0; i<s; i++)
      {
        Vector<number_type> zero(n,number_type(0));
        Z[i] = zero;
        for (int j = 0; j < n; j++)
        {
          Z[i][j] = zij[i*n+j];
        }
      }
      if (last_row_eq_b)
      {
        u += Z[s-1];
      }
      else
      {
        // compute ki
        Vector<number_type> zero(n,number_type(0));
        for (int i = 0; i < s; i++)
        {
          K[i] = zero;
          for (int j=0; j < s; j++)
          {
            K[i].update(Ainv[i][j],Z[j]);
          }
          K[i]*= (1.0/dt);

          // compute u
          u.update(dt*b[i], K[i]);
        }
      }
    }
      t = t+dt;
   }

   //! set current state
   void set_state (time_type t_, const Vector<number_type>& u_)
   {
     t = t_;
     u = u_;
   }

   //! get current state
   const Vector<number_type>& get_state () const
   {
     return u;
   }

   //! get current time
   time_type get_time () const
   {
     return t;
   }

   //! get dt used in last step (i.e. to compute current state)
   time_type get_dt () const
   {
     return dt;
   }

   //! how much should the ODE solver talk
   void set_verbosity(int verbosity_)
   {
     verbosity = verbosity_;
   }

  private:
    const M& model;
    time_type t, dt;
    Vector<number_type> u;
    Vector<number_type> w;
    Vector<Vector<number_type> > K;                     // save ki
    int n;											    // dimension of matrix A
    int s;
    DenseMatrix<number_type> A;				            // A, b, c as in the butcher tableau
	Vector<number_type> b;
	Vector<number_type> c;
    number_type sigma;
    int verbosity;
  };


  /** @brief Test convergence order of an ODE solver applied to a model problem

      \tparam M Type of model
      \tparam S Type of ODE solver

      \param model Model problem
      \param solver ODE solver
      \param T Solve to time T
      \param dt Roughest time step size
      \param l Number of different time step sizes dt, dt/2, dt/4, ...
   */
  template<class M, class S>
  void ordertest(const M& model,
                 S solver,
                 typename M::number_type T,
                 typename M::number_type h_0,
                 int l)
  {
    // Get types
    typedef typename M::time_type time_type;
    typedef typename M::number_type number_type;

    // error_array[i] = ||u(T)-u_i(T)||
    number_type error_array[l];

    Vector<number_type> exact_solution;
    model.exact_solution(T, exact_solution);

    for (int i=0; i<l; i++)
    {
      // Set initial time and value
      time_type t_start;
      Vector<number_type> initial_solution(1);
      model.initialize(t_start, initial_solution);
      solver.set_state(t_start, initial_solution);

      // Initial time step
      time_type dt = h_0/pow(2,i) ;
      solver.set_dt(dt);

      // Time loop
      while (solver.get_time()<T-2*solver.get_dt())
      {
        solver.step();
      }

      // Last steps
      if (solver.get_time()<T-solver.get_dt())
      {
        solver.set_dt((T-solver.get_time())/2.0);
        for(int i=0; i<2; i++)
        {
          solver.step();
        }
      }
      else
      {
        solver.set_dt(T-solver.get_time());
        solver.step();
      }

      // Error
      Vector<number_type> state = solver.get_state();
      error_array[i] = norm(exact_solution-state);

      if(i==0)
      {
        std::cout << "dt: "
                  << std::scientific << std::showpoint << std::setprecision(8)
                  << dt
                  << "  "
                  << "Error: "
                  << error_array[0] << std::endl;
      }
      if(i>0)
      {
        number_type rate = log(error_array[i-1]/error_array[i])/log(2);
        std::cout << "dt: "
                  << std::scientific << std::showpoint << std::setprecision(8)
                  << dt
                  << "  "
                  << "Error: "
                  << error_array[i]
                  << "  "
                  << "Rate: "
                  << rate << std::endl;
      }
    }
  }

} // namespace hdnum

#endif
