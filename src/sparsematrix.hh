// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/*
 * File:   sparsematrix.hh
 * Author: Christian Heusel <christian@heusel.eu>
 *
 * Created on August 25, 2020
 */

#ifndef SPARSEMATRIX_HH
#define SPARSEMATRIX_HH

#include <iomanip>
#include <iostream>
#include <vector>

#include "densematrix.hh"
#include "vector.hh"

namespace hdnum {

/*! \brief Sparse matrix Class with mathematical matrix operations
 */
template <typename REAL>
class SparseMatrix {
public:
    /** \brief Types used for array indices */
    using size_type = std::size_t;
    using VType = typename std::vector<REAL>;
    using VectorIterator = typename VType::iterator;
    using ConstVectorIterator = typename VType::const_iterator;

private:
    VType m_data;      // Matrix data is stored in an STL vector!
    size_type m_rows;  // Number of Matrix rows
    size_type m_cols;  // Number of Matrix columns

    static bool bScientific;
    static size_type nIndexWidth;
    static size_type nValueWidth;
    static size_type nValuePrecision;

    //! get matrix element for write access:
    REAL at(const size_type row, const size_type col) {}

    //! get matrix element for read-only access:
    const REAL at(const size_type row, const size_type col) const {
        return m_data[row * m_cols + col];
    }

public:
    //! default constructor (empty Matrix)
    SparseMatrix() : m_data(0, 0), m_rows(0), m_cols(0) {}

    //! constructor
    SparseMatrix(const size_type _rows, const size_type _cols,
                 const REAL def_val = 0)
        : m_data(_rows * _cols, def_val), m_rows(_rows), m_cols(_cols) {}

    //! constructor from initializer list
    SparseMatrix(const std::initializer_list<std::initializer_list<REAL>>& v) {}

    //! constructor from hdnum::DenseMatrix
    SparseMatrix(const hdnum::DenseMatrix<REAL>& other) {}

    void addNewRow(const hdnum::Vector<REAL>& rowvector) {}

    size_type rowsize() const { return m_rows; }
    size_type colsize() const { return m_cols; }

    // pretty-print output properties
    bool scientific() const { return bScientific; }

    // regular (possibly modifying) Iterators
    VectorIterator begin() {}
    VectorIterator end() {}

    // const Iterators
    ConstVectorIterator cbegin() const {}
    ConstVectorIterator cend() const {}
    ConstVectorIterator begin() const { return this->cbegin(); }
    ConstVectorIterator end() const { return this->cend(); }

    /*!
      \brief   Switch between floating point (default=true) and fixed point
      (false) display

      \b Example:
      \code
      hdnum::DenseMatrix<double> A(4,4);
      A.scientific(false); // fixed point representation for all DenseMatrix
      objects A.width(8); A.precision(3); identity(A);  // Defines the identity
      matrix of the same dimension std::cout << "A=" << A << std::endl; \endcode

      \b Output:
      \verbatim
      A=
      0        1        2        3
      0     1.000    0.000    0.000    0.000
      1     0.000    1.000    0.000    0.000
      2     0.000    0.000    1.000    0.000
      3     0.000    0.000    0.000    1.000
      \endverbatim
    */
    void scientific(bool b) const { bScientific = b; }

    //! get index field width for pretty-printing
    size_type iwidth() const { return nIndexWidth; }

    //! get data field width for pretty-printing
    size_type width() const { return nValueWidth; }

    //! get data precision for pretty-printing
    size_type precision() const { return nValuePrecision; }

    //! set index field width for pretty-printing
    void iwidth(size_type i) const { nIndexWidth = i; }

    //! set data field width for pretty-printing
    void width(size_type i) const { nValueWidth = i; }

    //! set data precision for pretty-printing
    void precision(size_type i) const { nValuePrecision = i; }

    // write access on matrix element A_ij using A(i,j)
    REAL& operator()(const size_type row, const size_type col) {}

    //! read-access on matrix element A_ij using A(i,j)
    const REAL& operator()(const size_type row, const size_type col) const {}

    //! read-access on matrix element A_ij using A[i][j]
    const ConstVectorIterator operator[](const size_type row) const {}

    //! write-access on matrix element A_ij using A[i][j]
    VectorIterator operator[](const size_type row) {}

    SparseMatrix operator=(const SparseMatrix& A) {}
    SparseMatrix operator=(const REAL value) {}

    bool operator==(const SparseMatrix& other) const {}
    bool operator!=(const SparseMatrix& other) const {}
    bool operator==(const hdnum::DenseMatrix<REAL>& other) const {}
    bool operator!=(const hdnum::DenseMatrix<REAL>& other) const {}

    // delete all the invalid comparisons
    bool operator<(const SparseMatrix& other) = delete;
    bool operator>(const SparseMatrix& other) = delete;
    bool operator<=(const SparseMatrix& other) = delete;
    bool operator>=(const SparseMatrix& other) = delete;

    SparseMatrix transpose() const {
        SparseMatrix A(m_cols, m_rows);
        for (size_type i = 0; i < m_rows; i++)
            for (size_type j = 0; j < m_cols; j++)
                A[j][i] = this->operator()(i, j);
        return A;
    }

    // Basic Matrix Operations
    [[nodiscard]] SparseMatrix operator+=(const SparseMatrix& B) {}
    [[nodiscard]] SparseMatrix operator-=(const SparseMatrix& B) {}
    [[nodiscard]] SparseMatrix operator*=(const REAL s) {}
    [[nodiscard]] SparseMatrix operator/=(const REAL s) {}

    void update(const REAL s, const SparseMatrix& B) {}

    template <class V>
    void mv(Vector<V>& y, const Vector<V>& x) const {}

    template <class V>
    void umv(Vector<V>& y, const Vector<V>& x) const {}

    template <class V>
    void umv(Vector<V>& y, const V& s, const Vector<V>& x) const {}

    void mm(const SparseMatrix<REAL>& A, const SparseMatrix<REAL>& B) {}

    [[nodiscard]] Vector<REAL> operator*(const Vector<REAL>& x) const {}

    [[nodiscard]] SparseMatrix operator+(const SparseMatrix& x) const {}
    [[nodiscard]] SparseMatrix operator-(const SparseMatrix& x) const {}
    [[nodiscard]] SparseMatrix operator*(const SparseMatrix& x) const {}
    [[nodiscard]] SparseMatrix operator/(const SparseMatrix& x) const {}

    //! compute row sum norm
    REAL norm_infty() const {
        REAL norm(0.0);
        for (size_type i = 0; i < rowsize(); i++) {
            REAL sum(0.0);
            for (size_type j = 0; j < colsize(); j++)
                sum += myabs((*this)(i, j));
            if (sum > norm) norm = sum;
        }
        return norm;
    }

    //! compute column sum norm
    REAL norm_1() const {
        REAL norm(0.0);
        for (size_type j = 0; j < colsize(); j++) {
            REAL sum(0.0);
            for (size_type i = 0; i < rowsize(); i++)
                sum += myabs((*this)(i, j));
            if (sum > norm) norm = sum;
        }
        return norm;
    }

    SparseMatrix<REAL> matchingIdentity() const {}
    static SparseMatrix identity(const size_type dimN) {}
};

template <typename REAL>
bool SparseMatrix<REAL>::bScientific = true;
template <typename REAL>
std::size_t SparseMatrix<REAL>::nIndexWidth = 10;
template <typename REAL>
std::size_t SparseMatrix<REAL>::nValueWidth = 10;
template <typename REAL>
std::size_t SparseMatrix<REAL>::nValuePrecision = 3;

template <typename REAL>
std::ostream& operator<<(std::ostream& out, const SparseMatrix<REAL>& A) {
    return out;
}

//! make a zero matrix
template <typename REAL>
inline void zero(SparseMatrix<REAL>& A) {}

/*!
  \relates SparseMatrix
  \n
  \b Function: make identity matrix
  \code
  template<class T>
  inline void identity (SparseMatrix<T> &A)
  \endcode
  \param[in] A reference to a SparseMatrix that shall be filled with entries

  \b Example:
  \code
  hdnum::SparseMatrix<double> A(4,4);
  identity(A);

  A.scientific(false); // fixed point representation for all DenseMatrix objects
  A.width(10);
  A.precision(5);

  std::cout << "A=" << A << std::endl;
  \endcode

  \b Output:
  \verbatim
  A=
  0          1          2          3
  0     1.00000    0.00000    0.00000    0.00000
  1     0.00000    1.00000    0.00000    0.00000
  2     0.00000    0.00000    1.00000    0.00000
  3     0.00000    0.00000    0.00000    1.00000
  \endverbatim

*/
template <class T>
inline void identity(SparseMatrix<T>& A) {}

}  // namespace hdnum

#endif  // SPARSEMATRIX_HH
