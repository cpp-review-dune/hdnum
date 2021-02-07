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
#include <map>
#include <numeric>
#include <string>
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
    using VectorIterator = typename std::vector<REAL>::iterator;
    using ConstVectorIterator = typename VType::const_iterator;

private:
    // Matrix data is stored in an STL vector!
    VType _data;

    // The non-null indices are stored in STL vectors with the size_type!
    // Explanation on how the mapping works can be found here:
    // https://de.wikipedia.org/wiki/Compressed_Row_Storage
    std::vector<size_type> _colIndices;
    std::vector<size_type> _rowPtr;

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
        return _data[row * m_cols + col];
    }

    // !function that converts container contents into
    // { 1, 2, 3, 4 }
    template <typename T>
    std::string comma_fold(T container) {
        return "{ " +
               std::accumulate(
                   std::next(container.begin()), container.end(),
                   std::to_string(container[0]),  // start with first element
                   [](std::string a, REAL b) {
                       return a + ", " + std::to_string(b);
                   }) +
               " }";
    };

public:
    //! default constructor (empty Matrix)
    SparseMatrix() noexcept
        : _data(), _colIndices(), _rowPtr(), m_rows(0), m_cols(0) {}

    //! constructor
    SparseMatrix(const size_type _rows, const size_type _cols)
        : _data(), _colIndices(), _rowPtr(_rows + 1), m_rows(_rows),
          m_cols(_cols) {}

    //! constructor from initializer list
    SparseMatrix(const std::initializer_list<std::initializer_list<REAL>> &v) {}

    size_type rowsize() const { return m_rows; }
    size_type colsize() const { return m_cols; }

    // pretty-print output properties
    bool scientific() const { return bScientific; }

    class column_iterator {
    public:
        using self_type = column_iterator;

        // conform to the iterator traits
        // https://en.cppreference.com/w/cpp/iterator/iterator_traits
        using difference_type = std::ptrdiff_t;
        using value_type = std::pair<REAL, size_type>;
        using pointer = value_type *;
        using reference = value_type &;
        using iterator_category = std::forward_iterator_tag;

        column_iterator(VectorIterator valIter,
                        std::vector<size_type>::iterator colIndicesIter)
            : _valIter(valIter), _colIndicesIter(colIndicesIter) {}

        // prefix
        self_type &operator++() {
            _valIter++;
            _colIndicesIter++;
            return *this;
        }

        // postfix
        self_type &operator++(int junk) {
            self_type cached = *this;
            _valIter++;
            _colIndicesIter++;
            return cached;
        }

        value_type operator*() {
            return std::make_pair(*_valIter, *_colIndicesIter);
        }
        value_type operator->() {
            return std::make_pair(*_valIter, *_colIndicesIter);
        }
        /* value_type operator*() { return std::make_pair(1, 2); } */
        /* value_type operator->() { return std::make_pair(3, 4); } */

        bool operator==(const self_type &other) {
            return (_valIter == other._valIter) and
                   (_colIndicesIter == other._colIndicesIter);
        }
        bool operator!=(const self_type &other) {
            return (_valIter != other._valIter) and
                   (_colIndicesIter != other._colIndicesIter);
        }

    private:
        VectorIterator _valIter;
        std::vector<size_type>::iterator _colIndicesIter;
    };

    class row_iterator {
    public:
        using self_type = row_iterator;

        // conform to the iterator traits
        // https://en.cppreference.com/w/cpp/iterator/iterator_traits
        using difference_type = std::ptrdiff_t;
        using value_type = VectorIterator;
        using pointer = VectorIterator *;
        using reference = VectorIterator &;
        using iterator_category = std::forward_iterator_tag;

        row_iterator(std::vector<size_type>::iterator rowPtrIter,
                     std::vector<size_type>::iterator colIndicesIter,
                     VectorIterator valIter)
            : _rowPtrIter(rowPtrIter), _colIndicesIter(colIndicesIter),
              _valIter(valIter) {}

        /* [[nodiscard]] VectorIterator begin() { return _valIter +
         * *_rowPtrIter; } */
        /* [[nodiscard]] VectorIterator end() { */
        /*     return _valIter + *(_rowPtrIter + 1); */
        /* } */

        [[nodiscard]] column_iterator begin() {
            return column_iterator((_valIter + *_rowPtrIter),
                                   (_colIndicesIter + *_rowPtrIter));
        }
        [[nodiscard]] column_iterator end() {
            return column_iterator((_valIter + *(_rowPtrIter + 1)),
                                   (_colIndicesIter + *(_rowPtrIter + 1)));
        }

        // prefix
        self_type &operator++() {
            _rowPtrIter++;
            _currRow++;
            return *this;
        }

        // postfix
        self_type &operator++(int junk) {
            self_type cached = *this;
            _rowPtrIter++;
            _currRow++;
            return cached;
        }

        self_type operator*() { return *this; }
        self_type operator->() { return *this; }
        bool operator==(const self_type &rhs) {
            return _rowPtrIter == rhs._rowPtrIter;
        }
        bool operator!=(const self_type &rhs) {
            return _rowPtrIter != rhs._rowPtrIter;
        }

    private:
        std::vector<size_type>::iterator _rowPtrIter;
        size_type _currRow;
        std::vector<size_type>::iterator _colIndicesIter;
        VectorIterator _valIter;
    };

    /* class const_iterator { */
    /* public: */
    /*     using self_type = const_iterator; */

    /*     // conform to the iterator traits */
    /*     // https://en.cppreference.com/w/cpp/iterator/iterator_traits */
    /*     using difference_type = std::ptrdiff_t; */
    /*     using value_type = REAL; */
    /*     using pointer = REAL *; */
    /*     using reference = REAL &; */
    /*     using iterator_category = std::forward_iterator_tag; */

    /*     const_iterator(pointer ptr) : _ptr {ptr} {} */

    /*     // prefix */
    /*     self_type operator++() { */
    /*         _ptr++; */
    /*         return *this; */
    /*     } */

    /*     // postfix */
    /*     self_type operator++(int junk) { */
    /*         self_type cached = *this; */
    /*         _ptr++; */
    /*         return cached; */
    /*     } */

    /*     const value_type &operator*() { return *_ptr; } */
    /*     const pointer operator->() { return _ptr; } */
    /*     bool operator==(const self_type &rhs) { return _ptr == rhs._ptr; } */
    /*     bool operator!=(const self_type &rhs) { return _ptr != rhs._ptr; } */

    /* private: */
    /*     pointer _ptr; */
    /* }; */

    // regular (possibly modifying) Iterators
    [[nodiscard]] row_iterator begin() {
        return SparseMatrix<REAL>::row_iterator(
            _rowPtr.begin(), _colIndices.begin(), _data.begin());
    }
    [[nodiscard]] row_iterator end() {
        return row_iterator(_rowPtr.end() - 1, _colIndices.begin(),
                            _data.begin());
    }

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
    REAL &operator()(const size_type row, const size_type col) {
        _data.push_back(REAL {});
        _colIndices.push_back(col);
        if (row < _rowPtr[col]) {
            _rowPtr[col] = row;
        }
        return _data[_data.size() - 1];
    }

    //! read-access on matrix element A_ij using A(i,j)
    const REAL &operator()(const size_type row, const size_type col) const {
        if (m_cols - 1 <= col) {
            HDNUM_ERROR("Out of bounds access: column too big!");
        } else if (m_rows - 1 <= row) {
            HDNUM_ERROR("Out of bounds access: row too big!");
        }

        // Handle the zero matrix case
        if (_rowPtr[col] == 0 and _rowPtr[col + 1] == 0 or
            _rowPtr.size() == 1) {
            return REAL {};
        }

        // look for the entry
        for (auto i = _rowPtr[col]; i < _rowPtr[col + 1]; ++i) {
            if (_colIndices[i] == col) {
                return _data[i];
            }
        }
        // look for the entry
        return REAL {};
    }

    //! read-access on matrix element A_ij using A[i][j]
    const ConstVectorIterator operator[](const size_type row) const {}

    //! write-access on matrix element A_ij using A[i][j]
    VectorIterator operator[](const size_type row) {}

    SparseMatrix &operator=(const SparseMatrix &other) {
        _data = other._data;
        _rowPtr = other._rowPtr;
        _colIndices = other._colIndices;

        m_cols = other.m_cols;
        m_rows = other.m_rows;
        return *this;
    }

    [[nodiscard]] bool operator==(const SparseMatrix &other) const {
        return (_data == other._data) and              //
               (_rowPtr == other._rowPtr) and          //
               (_colIndices == other._colIndices) and  //
               (m_cols == other.m_cols) and            //
               (m_rows == other.m_rows);
    }
    [[nodiscard]] bool operator!=(const SparseMatrix &other) const {
        return not (*this == other);
    }

    // delete all the invalid comparisons
    bool operator<(const SparseMatrix &other) = delete;
    bool operator>(const SparseMatrix &other) = delete;
    bool operator<=(const SparseMatrix &other) = delete;
    bool operator>=(const SparseMatrix &other) = delete;

    SparseMatrix transpose() const {
        SparseMatrix A(m_cols, m_rows);
        for (size_type i = 0; i < m_rows; i++)
            for (size_type j = 0; j < m_cols; j++)
                A[j][i] = this->operator()(i, j);
        return A;
    }

    // Basic Matrix Operations
    [[nodiscard]] SparseMatrix operator+=(const SparseMatrix &B) {}
    [[nodiscard]] SparseMatrix operator-=(const SparseMatrix &B) {}
    [[nodiscard]] SparseMatrix operator*=(const REAL s) {}
    [[nodiscard]] SparseMatrix operator/=(const REAL s) {}

    void update(const REAL s, const SparseMatrix &B) {}

    template <class V>
    void mv(Vector<V> &y, const Vector<V> &x) const {}

    template <class V>
    void umv(Vector<V> &y, const Vector<V> &x) const {}

    template <class V>
    void umv(Vector<V> &y, const V &s, const Vector<V> &x) const {}

    void mm(const SparseMatrix<REAL> &A, const SparseMatrix<REAL> &B) {}

    [[nodiscard]] Vector<REAL> operator*(const Vector<REAL> &x) const {}

    [[nodiscard]] SparseMatrix operator+(const SparseMatrix &x) const {}
    [[nodiscard]] SparseMatrix operator-(const SparseMatrix &x) const {}
    [[nodiscard]] SparseMatrix operator*(const SparseMatrix &x) const {}
    [[nodiscard]] SparseMatrix operator/(const SparseMatrix &x) const {}

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

    std::string to_string() noexcept {
        return "values=" + comma_fold(_data) + "\n" +        //
               "colInd=" + comma_fold(_colIndices) + "\n" +  //
               "rowPtr=" + comma_fold(_rowPtr) + "\n";       //
    }

    void print() noexcept { std::cout << this->to_string(); }

    SparseMatrix<REAL> matchingIdentity() const {}
    static SparseMatrix identity(const size_type dimN) {}

    class builder {
        size_type m_rows {};  // Number of Matrix rows, 0 by default
        size_type m_cols {};  // Number of Matrix columns, 0 by default
        std::vector<std::map<size_type, REAL>> _rows;

    public:
        builder(size_type new_m_rows, size_type new_m_cols)
            : m_rows {new_m_rows}, m_cols {new_m_cols}, _rows {m_rows} {}

        builder() = default;

        std::pair<typename std::map<size_type, REAL>::iterator, bool> addEntry(
            size_type i, size_type j, REAL value) {
            return _rows.at(i).emplace(j, value);
        }

        std::pair<typename std::map<size_type, REAL>::iterator, bool> addEntry(
            size_type i, size_type j) {
            return addEntry(i, j, REAL {});
        };

        [[nodiscard]] bool operator==(
            const SparseMatrix::builder &other) const {
            return (m_rows == other.m_rows) and (m_cols == other.m_cols) and
                   (_rows == other._rows);
        }

        [[nodiscard]] bool operator!=(
            const SparseMatrix::builder &other) const {
            return not (*this == other);
        }

        [[nodiscard]] size_type colsize() noexcept { return m_cols; }
        [[nodiscard]] size_type rowsize() noexcept { return m_rows; }

        size_type setNumCols(size_type new_m_cols) noexcept {
            m_cols = new_m_cols;
            return m_cols;
        }
        size_type setNumRows(size_type new_m_rows) {
            m_rows = new_m_rows;
            _rows.resize(m_cols);
            return m_rows;
        }

        void clear() noexcept {
            for (auto &row : _rows) {
                row.clear();
            }
        }

        [[nodiscard]] std::string to_string() {
            std::string output;
            for (std::size_t i = 0; i < _rows.size(); i++) {
                for (const auto &[index, value] : _rows[i]) {
                    output += "i=" + std::to_string(i) +
                              ", j=" + std::to_string(index) + " => " +
                              std::to_string(value) + "\n";
                }
            }
            return output;
        }

        [[nodiscard]] SparseMatrix build() {
            auto result = SparseMatrix<REAL>(m_rows, m_cols);

            for (std::size_t i = 0; i < _rows.size(); i++) {
                result._rowPtr[i + 1] = result._rowPtr[i];
                for (const auto &[index, value] : _rows[i]) {
                    result._colIndices.push_back(index);
                    result._data.push_back(value);
                    result._rowPtr[i + 1]++;
                }
            }
            return result;
        }
    };
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
std::ostream &operator<<(std::ostream &out, const SparseMatrix<REAL> &A) {
    return out;
}

//! make a zero matrix
template <typename REAL>
inline void zero(SparseMatrix<REAL> &A) {}

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
  // fixed point representation for all DenseMatrix objects
  A.scientific(false);
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
inline void identity(SparseMatrix<T> &A) {}

}  // namespace hdnum

#endif  // SPARSEMATRIX_HH
