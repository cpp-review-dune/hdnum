#ifndef HDNUM_MODELPROBLEM_HIGH_DIM_HH
#define HDNUM_MODELPROBLEM_HIGH_DIM_HH
/** @brief Example class for a differential equation model

    The model is

		u'(t) = (5 -2// -2 5)*u(t)
		u(0) = (1 3)
    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/

template<class T, class N=T>
class ModelProblem_high_dim
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor
  ModelProblem_high_dim ()
  {}

  //! return number of components for the model
  std::size_t size () const
  {
    return 2; // as the dimension of the differential equation ist two
  }

  //! set initial state including time value
  void initialize (T& t0, hdnum::Vector<N>& x0) const
  {
    t0 = 0;
    x0[0] = 1.0;
		x0[1] = 3.0;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = 5*x[0] - 2*x[1];
		result[1] = -2*x[0]+5*x[1];
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
		result[0][0] = 5.0;
		result[0][1] = -2.0;
		result[1][0] = -2.0;
		result[1][1] = 5.0;
  }

private:
  //
};
#endif
