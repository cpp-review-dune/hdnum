#include <iostream>
#include <vector>
#include "hdnum.hh"


template<class N>
class WurzelProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  WurzelProblem (number_type a_)
    : a(a_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 1;
  }

  //! model evaluation
  void F (const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = x[0]*x[0] - a;
  }

  //! jacobian evaluation needed for implicit solvers
  void F_x (const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
    result[0][0] = number_type(2.0)*x[0];
  }

private:
  number_type a;
};



int main ()
{
  //typedef double Number;                 // numbertype
  typedef mpf_class Number;
  mpf_set_default_prec(1024);

  typedef WurzelProblem<Number> Problem; // problem type
  Problem problem(2.0);                  // an instance of the problem

  hdnum::Banach banach;                         // an Object of class Banach
  banach.set_maxit(20);                  // set parameters
  banach.set_verbosity(2);
  banach.set_reduction(1e-100);
  banach.set_abslimit(1e-100);
  banach.set_linesearchsteps(3);
  banach.set_sigma(0.1);

  hdnum::Vector<Number> u(problem.size());      // vector for storing solvations
  u[0] = 2.0;                            // value for starting
  banach.solve(problem,u);               // compute solvation
  std::cout << "Ergebnis: " << u[0] << std::endl;

  return 0;
}
