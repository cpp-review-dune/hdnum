// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_TESTODE_HH
#define HDNUM_TESTODE_HH

#include<vector>
#include "../../src/newton.hh"

/** @file
 *  @brief solvers for ordinary differential equations
 */

namespace hdnum {



/** @brief classical Runge-Kutta method (order n with n stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class RungeKutta_n
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model - ich gehe davon aus, dass Dimension von Matrix und Vektor, sowie der Zahlentyp zusammenpasst
    RungeKutta_n (const M& model_, DenseMatrix<double> Mat, Vector<double> BV, Vector<double> CV)
      : model(model_), u(model.size()), w(model.size())
    {
      A = Mat;
      B = BV;
      C = CV;
      n = Mat.rowsize ();
      model.initialize(t,u);
      dt = 0.1;
      K [n];              // ein Array der Größe n erzeugen
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {

      // k1 berechnen
        w = u;
        model.f(t, w, K[0]);

        for (int i = 0; i < n; i++)
        {
            number_type sum = K[0]*B[0];
            for (int j = 1; j < i-1; j++)       //berechne ki
            {
                sum = sum + A[i-1][j-1]*K[j-1];
            }
            model.f(t + C[i-1]*dt, w + dt*sum, K[i-1]);
            u.update(dt *B[i-1], K[i-1]);
        }



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
    const M& model;
    time_type t, dt;
    Vector<number_type> u,w;
    Vector<number_type> K;                          // save ki
    int n;											// dimension of matrix A
    DenseMatrix<double> A;				            // A, B, C as in the butcher tableau
	Vector<double> B;
	Vector<double> C;
  };

} // namespace hdnum

#endif
