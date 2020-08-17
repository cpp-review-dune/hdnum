#include <iostream>
#include <vector>
#include "hdnum.hh"


template<class N>
class EllipseProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  EllipseProblem ()
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 2;
  }

  //! model evaluation
  void F (const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = x[0]*x[0] + x[1]*x[1] - 4.0;
    result[1] = x[0]*x[0]/9.0 + x[1]*x[1] - 1.0;
  }

  //! jacobian evaluation needed for implicit solvers
  void F_x (const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
    result[0][0] = number_type(2.0)*x[0]; result[0][1] = number_type(2.0)*x[1];
    result[1][0] = number_type(2.0)*x[0]/number_type(9.0); result[1][1] = number_type(2.0)*x[1];
  }
};



int main ()
{
  typedef double Number;                 // numbertype
  // typedef mpf_class Number;
  // mpf_set_default_prec(1024);

  typedef EllipseProblem<Number> Problem; // problem type
  Problem problem;                  // an instance of the problem

  hdnum::Banach banach;                         // an Object of class Banach
  banach.set_maxit(5000);                  // set parameters
  banach.set_verbosity(2);
  banach.set_reduction(1e-10);
  banach.set_abslimit(1e-20);
  banach.set_sigma(-0.01);

  hdnum::Vector<Number> u(problem.size());   // vector for storing solvations
  u[0] = -2.0; u[1] = 1.0;                    // value for starting
  // u[0] = -1.83711730; u[1] = 0.7905694151;                    // value for starting

  hdnum::Newton newton;                         // Ein Newtonobjekt
  newton.set_maxit(20);                  // Setze diverse Parameter
  newton.set_verbosity(2);
  newton.set_reduction(1e-10);
  newton.set_abslimit(1e-20);
  newton.set_linesearchsteps(3);


  //banach.solve(problem,u);               // compute solvation
  newton.solve(problem,u);               // compute solvation
  std::cout << "Ergebnis: " <<  std::setprecision(15) << u[0] << " " << u[1] << std::endl;


  return 0;
}
