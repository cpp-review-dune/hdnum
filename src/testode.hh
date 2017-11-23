// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_TESTODE_HH
#define HDNUM_TESTODE_HH

#include<vector>
#include "newton.hh"

/** @file
 *  @brief solvers for ordinary differential equations
 */

namespace hdnum {

template<class T>
class compute_K
{
public:
    /** \brief export size_type */
    typedef typename T::size_type size_type;

    /** \brief export time_type */
    typedef typename T::time_type time_type;

    /** \brief export number_type */
    typedef typename T::number_type number_type; 

  //! constructor stores parameter lambda
  compute_K (const T& model_, DenseMatrix<double> Mat, Vector<double> BV, Vector<double> CV)
      : model(model_), u(model.size()), w(model.size()), K(Mat.rowsize ())
    {
      A = Mat;
      B = BV;
      C = CV;
      s = Mat.rowsize ();
      model.initialize(t,u);
      dt = 0.1;
      n = model.size();
      for (int i = 0; i<s; i++){
        K[i].resize(n);
      }
    }

  //! return number of componentes for the model
  std::size_t size () const
  {
    return n*s;
  }

  //! model evaluation
  void F (const Vector<number_type>& x, Vector<number_type>& result) const
  {
    Vector<Vector<number_type>> xx (s);    //hilfsvektor
    
    for (int i = 0; i < s; i++)
    {
        xx[i].resize(n);
        for(int k = 0; k < n; k++)
        {
            xx[i][k] = x[i*n + k];
        }
    }
    
    Vector<Vector<number_type>> f (s);
    for (int i = 0; i < s; i++)
    {
        f[i].resize(n);
        model.f(t + C[i] * dt, u + xx[i], f[i]);
    }


    Vector<Vector<number_type>> hr (s);
    for (int i = 0; i < s; i++)
    {
        hr[i].resize(n);
    }

    for (int i = 0; i < s; i++)
    {
        Vector<number_type> sum (n, 0.0);
        
        for (int j = 0; j < s; j++)
        {
            sum.update(dt*A[i][j], f[j]);
        }
        

        hr[i]  = xx[i] - sum;        

    }

    //uebersetze hr nach result
    for (int i = 0; i< s; i++)
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
    Vector<Vector<number_type>> xx (s);    //hilfsvektor
    
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
            DenseMatrix<number_type> J (n, n, 0.0);
            DenseMatrix<number_type> H (n, n, 0.0);     //Hilfsmatrix
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
    Vector<number_type> w;
    Vector<Vector<number_type>> K;                          // save ki
    int n, s;											// dimension of matrix A
    DenseMatrix<double> A;				            // A, B, C as in the butcher tableau
	Vector<double> B;
	Vector<double> C;
};


/** @brief classical Runge-Kutta method (order n with n stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class N>
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
    RungeKutta_n (const N& model_, DenseMatrix<double> Mat, Vector<double> BV, Vector<double> CV)
      : model(model_), u(model.size()), w(model.size()), K(Mat.rowsize ())
    {
      A = Mat;
      B = BV;
      C = CV;
      n = Mat.rowsize ();
      model.initialize(t,u);
      dt = 0.1;
      for (int i = 0; i<n; i++){
        K[i].resize(model.size());
      }
      //K [n];              // ein Array der Größe n erzeugen
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    bool check_explicit ()
    {
        bool ergebnis = true;
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
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
          // k1 berechnen
            w = u;

            model.f(t, w, K[0]);        

            for (int i = 0; i < n; i++)
            {
                Vector<number_type> sum (K[0].size(), 0.0);
                sum.update(B[0], K[0]);
                for (int j = 0; j < i+1; j++)       //berechne ki
                {
                    //sum = sum + A[i-1][j-1]*K[j-1];
                    sum.update(A[i][j],K[j]);
                }
                Vector<number_type> wert = w.update(dt,sum);

                model.f(t + C[i]*dt, wert, K[i]);

                u.update(dt *B[i], K[i]);
            }
        }
        else
        {
              compute_K<N> problem(model, A, B, C); // Problemtyp


              Banach banach;                         // Ein Banachobjekt
              banach.set_maxit(20);                  // Setze diverse Parameter
              banach.set_verbosity(0);
              banach.set_reduction(1e-100);
              banach.set_abslimit(1e-100);
              banach.set_linesearchsteps(3);
              banach.set_sigma(0.1);
              Vector<number_type> kij (model.size()*n,0.0);
              banach.solve(problem,kij);               // Berechne Lösung
              Vector<Vector<number_type>> K (n);
            for(int i = 0; i< n; i++)
            {
                K[i].resize(model.size());
                for (int j = 0; j< model.size(); j++)
                {
                    K[i][j] = kij[i*model.size()+j];
                }

                u.update(dt*B[i], K[i]);
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



  private:
    const N& model;
    time_type t, dt;
    Vector<number_type> u;
    Vector<number_type> w;
    Vector<Vector<number_type>> K;                          // save ki
    int n;											// dimension of matrix A
    DenseMatrix<double> A;				            // A, B, C as in the butcher tableau
	Vector<double> B;
	Vector<double> C;
  };

} // namespace hdnum

#endif
