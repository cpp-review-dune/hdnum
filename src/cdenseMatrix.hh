#ifndef CDenseMatrix_HH
#define CDenseMatrix_HH

#include "densematrix.hh"
#include "cvector.hh"

namespace hdnum {

  /*! \brief Class with mathematical matrix operations (for complex values)
   */
  template<typename N>
  class CDenseMatrix: public DenseMatrix<std::complex<N>>
  {
  public:
	/** \brief Type used for array indices */
	typedef std::size_t size_type;
	typedef typename std::vector<std::complex<N>> VType;
	typedef typename VType::const_iterator ConstVectorIterator;
	typedef typename VType::iterator VectorIterator;

  private:
  //! function to calculate the modulus of a value
  N myabs (std::complex<N> x) const
  {
    return sqrt(x.real() * x.real() + x.imag() * x.imag());
  }

  public:
	//! default constructor (empty Matrix)
	CDenseMatrix()
  :DenseMatrix<std::complex<N>>()
	{
	}

	//! constructor
	CDenseMatrix( const std::size_t _rows,
				 const std::size_t _cols,
				 const N def_val=0
				 )
         :DenseMatrix<std::complex<N>>(_rows,_cols, def_val)
	{
	}

  //! constructor from initializer list
  CDenseMatrix (const std::initializer_list<std::initializer_list<std::complex<N>>> &v)
  :DenseMatrix<std::complex<N>>(v)
  {
  }

  //! constructor from complex Densematrix
  CDenseMatrix(DenseMatrix<std::complex<N>> A)
  {
    this->m_rows = A.rowsize();
    this->m_cols = A.colsize();
    for(int i = 0; i< A.rowsize();++i)
      for(int j = 0; j < A.colsize();++j)
        this->m_data.push_back(A[i][j]);
  }
    
  /*!
    \brief Conjugate Transposition

    Return the conjugate transposed as a new matrix.
  */
  CDenseMatrix ctranspose () const
  {
    CDenseMatrix A(this->m_cols,this->m_rows);
    const CDenseMatrix &self = *this;
    for (size_type i=0; i<this->m_rows; i++)
      for (size_type j=0; j<this->m_cols; j++)
        A[j][i] = std::conj(self[i][j]);
    return A;
  }

  //! compute row sum norm
  N norm_infty () const
  {
    N norm(0.0);
    for (std::size_t i=0; i<this->rowsize(); i++)
      {
        N sum(0.0);
        for (std::size_t j=0; j<this->colsize(); j++)
          sum += myabs((*this)(i,j));
        if (sum>norm) norm = sum;
      }
    return norm;
  }

  //! compute column sum norm
  N norm_1 () const
  {
    N norm(0.0);
    for (std::size_t j=0; j<this->colsize(); j++)
      {
        N sum(0.0);
        for (std::size_t i=0; i<this->rowsize(); i++)
          sum += myabs((*this)(i,j));
        if (sum>norm) norm = sum;
      }
    return norm;
  }

  };





}  // namespace hdnum

#endif  // CDenseMatrix_HH
