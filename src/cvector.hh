#ifndef CVector_HH
#define CVector_HH

#include "vector.hh"

namespace hdnum {
  /*! \brief Class with mathematical vector operations
   */

  // Vector for complex values
  template<typename N>
  class CVector : public Vector<std::complex<N>>  // inherit from the STL vector
  {
  public:
    /** \brief Type used for array indices */
    typedef std::size_t size_type;
    typedef N value_type;

  public:

    //! default constructor, also inherited from the STL vector default constructor
    CVector() : Vector<std::complex<N>>()
    {
    }

    //! another constructor, with arguments, setting the default value for all entries of the vector of given size
    CVector( const size_t size,                 // user must specify the size
            const std::complex<N> defaultvalue_ = {0,0}    // if not specified, the value 0 will take effect
            )
      : Vector<std::complex<N>>( size, defaultvalue_ )
    {
    }

    //! constructor from initializer list
    CVector (const std::initializer_list<N> &v)
    {
      for (auto elem : v) this->push_back(elem);
    }

    //! constructor from complex vector
    CVector(Vector<std::complex<N>> v){
      for(auto vec:v){
        this->push_back(vec);
      }
    }

    //! Square of the Euclidean norm
    N two_norm_2() const
    {
      N sum( 0 );
      const CVector & self = *this;
      for (size_t i = 0; i < (size_t) this->size(); ++i)
        sum += (std::conj(self[i]) * self[i]).real();
      return sum;
    }

    N two_norm() const
    {
      return sqrt(two_norm_2());
    }

  };

  //! norm of a complex vector
  template<class N>
  inline N norm (CVector<N> x)
  {
    N sum(0.0);
    for (typename CVector<N>::size_type i=0; i<x.size(); i++)
      sum += (std::conj(x[i])*x[i]).real();
    return sqrt(sum);
  }

  // Check if a vector is complex or not 
  struct ComplexChecker{
    template< typename N>
      constexpr bool isComplexVector(const Vector<N>& n) {
      return false;
    }
    template< typename N>
      constexpr bool isComplexVector(const CVector<N>& n) {
      return true;
    }
  };

} // end of namespace hdnum

#endif	/* _CVector_HH */
