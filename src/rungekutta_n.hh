// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_RUNGEKUTTA_N_HH
#define HDNUM_RUNGEKUTTA_N_HH

#include<vector>
#include "newton.hh"

/** @file
 *  @general Runge-Kutta solver
 */

  namespace hdnum {

  template<class T>
  class compute_Z
  {
  public:
    /** \brief export size_type */
    typedef typename T::size_type size_type;

    /** \brief export time_type */
    typedef typename T::time_type time_type;

    /** \brief export number_type */
    typedef typename T::number_type number_type; 

    //! constructor stores parameter lambda
    compute_Z (const T& model_, DenseMatrix<number_type> Mat, Vector<number_type> BV, Vector<number_type> CV, time_type t_, Vector<number_type> u_, time_type dt_)
        : model(model_) , u(model.size())
      {
        A = Mat;
        B = BV;
        C = CV;
        s = Mat.rowsize ();
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
      Vector<Vector<number_type>> xx (s);    
      for (int i = 0; i < s; i++)
      {
        xx[i].resize(n,number_type(0));
        for(int k = 0; k < n; k++)
        {
          xx[i][k] = x[i*n + k];
        }
      }
      Vector<Vector<number_type>> f (s);
      for (int i = 0; i < s; i++)
      {
        f[i].resize(n, number_type(0));
        model.f(t + C[i] * dt, u + xx[i], f[i]);
      }
      Vector<Vector<number_type>> hr (s);
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
      Vector<Vector<number_type>> xx (s);   
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
          model.f_x(t+C[j]*dt, u + xx[j],H);
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
    const T& model;
    time_type t, dt;
    Vector<number_type> u;
    int n, s;						    // dimension of matrix A and model.size
    DenseMatrix<number_type> A;				// A, B, C as in the butcher tableau
    Vector<number_type> B;
    Vector<number_type> C;
  };


/** @brief classical Runge-Kutta method (order n with n stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
*/
  template<class N, class S = Newton>
  class RungeKutta_n
  {
  public:
    /** \brief export size_type */
    typedef typename N::size_type size_type;

    /** \brief export time_type */
    typedef typename N::time_type time_type;

    /** \brief export number_type */
    typedef typename N::number_type number_type;

    //! constructor stores reference to the model - ich gehe davon aus, dass Dimension von Matrix und Vektor, sowie der Zahlentyp zusammenpasst
    RungeKutta_n (const N& model_, DenseMatrix<number_type> Mat, Vector<number_type> BV, Vector<number_type> CV)
      : model(model_), u(model.size()), w(model.size()), K(Mat.rowsize ())
    {
      A = Mat;
      B = BV;
      C = CV;
      s = Mat.rowsize ();
      n = model.size();
      model.initialize(t,u);
      dt = 0.1;
      for (int i = 0; i < s; i++)
      {
        K[i].resize(n, number_type(0));
      }
      sigma = 0.01;
      verbosity = 0;
   }

   //! constructor stores reference to the model - ich gehe davon aus, dass Dimension von Matrix und Vektor, sowie der Zahlentyp zusammenpasst, hier kann man auserdem noch Sigma setzen
   RungeKutta_n (const N& model_, DenseMatrix<number_type> Mat, Vector<number_type> BV, Vector<number_type> CV, number_type sigma_)
     : model(model_), u(model.size()), w(model.size()), K(Mat.rowsize ())
   {
     A = Mat;
     B = BV;
     C = CV;
     s = Mat.rowsize ();
     n = model.size();
     model.initialize(t,u);
     dt = 0.1;
     for (int i = 0; i < s; i++)
     {
       K[i].resize(n, number_type(0));
     }
     sigma = sigma_;
     verbosity = 0;
  }

  //! set time step for subsequent steps
  void set_dt (time_type dt_)
  {
    dt = dt_;
  }
  bool check_explicit ()
  {
    bool ergebnis = true;
    for (int i = 0; i < s; i++)
    {
      for (int j = i; j < s; j++)
      {
        if (A[i][j] != 0.0)
        {
          ergebnis = false;
        }
      }
    }
    return ergebnis;
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
        sum.update(B[0], K[0]);
        //compute k_i
        for (int j = 0; j < i+1; j++)       
        {
          sum.update(A[i][j],K[j]);
        }
        Vector<number_type> wert = w.update(dt,sum);
        model.f(t + C[i]*dt, wert, K[i]);
        u.update(dt *B[i], K[i]);
      }
    }
    if (not check_explicit()) 
    {
      compute_Z<N> problem(model, A, B, C, t, u, dt);           // problemtype
      bool last_row_eq_b = true;
      for (int i = 0; i<s; i++)
      {
        if (A[s-1][i] != B[i])     
        {
          last_row_eq_b = false;
        }
      }
      S Solver;                         
      Solver.set_maxit(2000);
      Solver.set_verbosity(verbosity);
      Solver.set_reduction(1e-10);
      Solver.set_abslimit(1e-10);
      Solver.set_linesearchsteps(10);
      Solver.set_sigma(0.01);
      Vector<number_type> zij (s*n,0.0);
      Solver.solve(problem,zij);                                // compute solution
      Vector<Vector<number_type>> Z (s, 0.0);
      DenseMatrix<number_type> Ainv (s,s,number_type(0));
      if (not last_row_eq_b) 
      // compute invers of A with LR decomposition
      {
        for (int i=0; i < s; i++)                       
        {
          Vector<number_type> e (s, number_type(0));
          e[i]=number_type(1);
          Vector<number_type> w (s, number_type(0));
          Vector<number_type> x (s, number_type(0));
          Vector<number_type> z (s, number_type(0));
          Vector<std::size_t> p(s);
          Vector<std::size_t> q(s);
          DenseMatrix<number_type> Temp (s,s,0.0);
          Temp = A;
          row_equilibrate(A,w);                         
          lr_fullpivot(A,p,q);                          
          apply_equilibrate(w,e);                       
          permute_forward(p,e);                         
          solveL(A,e,e);                  
          solveR(A,z,e);                                
          permute_backward(q,z);                        
          for (int j = 0; j < s; j++)
          {
	        Ainv[j][i] = z[j];
          }

          A = Temp;
        }
      }
      for(int i = 0; i < s; i++)
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
          u.update(dt*B[i], K[i]);      // wenn man dt hier weg laesst, sieht das ergebnis besser aus :(
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

   void set_verbosity(int verbosity_)
   {
     verbosity = verbosity_;
   }

  private:
    const N& model;
    time_type t, dt;
    Vector<number_type> u;
    Vector<number_type> w;
    Vector<Vector<number_type>> K;                      // save ki
    int n;											    // dimension of matrix A
    int s;
    DenseMatrix<number_type> A;				            // A, B, C as in the butcher tableau
	Vector<number_type> B;
	Vector<number_type> C;
    number_type sigma;
    int verbosity;
  };

  typedef double Number;
  template<class N, class S = Newton>
  void ordertest(const N& model, DenseMatrix<Number> Mat, Vector<Number> BV, Vector<Number> CV, Number T, Number h_0, int L)
  {
    /** \brief export size_type */
    typedef typename N::size_type size_type;

    /** \brief export time_type */
    typedef typename N::time_type time_type;

    /** \brief export number_type */
    typedef typename N::number_type number_type;

    // aim U[i] = u_i(T)
    Vector<number_type> U[L];
    // aim Earray[i] = ||u(T)-u_i(T)||
    number_type Earray[L];
    number_type alpha[L-1];
    Vector<number_type> exact_solution;
    model.u(T, exact_solution);
    for (int i = 0; i<L; i++)
    {
      RungeKutta_n<N,S> solver(model, Mat, BV, CV);
      solver.set_dt(h_0/pow(2,i));                      // set initial time step
      Vector<time_type> times;                          // store time values here
      Vector<Vector<number_type>> states;               // store states here
      times.push_back(solver.get_time());               // initial time
      states.push_back(solver.get_state());             // initial state
      while (solver.get_time()<T-2*solver.get_dt())     // the time loop
      {
        solver.step();                                  // advance model by one time step
        times.push_back(solver.get_time());             // save time
        states.push_back(solver.get_state());           // and state
      }
      // three different cases to get to end time T
      time_type t = times[times.size()-1];
      number_type eps = 0.000001;
      if (t < T-solver.get_dt())
      {
        solver.set_dt((T-t)/2.0);
        for(int i = 0; i < 2;i++)
        {
          solver.step();                                // advance model by one time step
          times.push_back(solver.get_time());           // save time
          states.push_back(solver.get_state());         // and state
        }
      }
      else if ( t < T and abs(t-T) > eps)
      {
        solver.set_dt(T-t);
        solver.step();                                  // advance model by one time step
        times.push_back(solver.get_time());             // save time
        states.push_back(solver.get_state());           // and state
      }
      Earray[i] = norm(exact_solution-states[states.size()-1]);
      if(i == 0)
      {
        std::cout << std::scientific << std::showpoint
                  << std::setprecision(8) << Earray[0] << std::endl;
      }
      if(i > 0)
      {
        alpha[i-1] = log(Earray[i-1]/Earray[i])/log(2);
        std::cout << std::scientific << std::showpoint
                  << std::setprecision(8) << Earray[i] << " " << alpha[i-1] << std::endl;
      }
    }
  }
} // namespace hdnum

#endif
