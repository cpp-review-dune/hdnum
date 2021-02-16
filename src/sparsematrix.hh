// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/*
 * File:   sparsematrix.hh
 * Author: Christian Heusel <christian@heusel.eu>
 *
 * Created on August 25, 2020
 */

#ifndef SPARSEMATRIX_HH
#define SPARSEMATRIX_HH

#include <algorithm>
#include <functional>
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
    static const REAL _zero;

    // !function that converts container contents into
    // { 1, 2, 3, 4 }
    template <typename T>
    [[nodiscard]] std::string comma_fold(T container) const {
        return "{ " +
               std::accumulate(
                   std::next(container.cbegin()), container.cend(),
                   std::to_string(container[0]),  // start with first element
                   [](std::string a, REAL b) {
                       return a + ", " + std::to_string(b);
                   }) +
               " }";
    }

public:
    //! default constructor (empty Matrix)
    SparseMatrix() noexcept
        : _data(), _colIndices(), _rowPtr(), m_rows(0), m_cols(0) {}

    //! constructor
    SparseMatrix(const size_type _rows, const size_type _cols)
        : _data(), _colIndices(), _rowPtr(_rows + 1), m_rows(_rows),
          m_cols(_cols) {}

    [[nodiscard]] size_type rowsize() const { return m_rows; }
    [[nodiscard]] size_type colsize() const { return m_cols; }

    // pretty-print output properties
    [[nodiscard]] bool scientific() const { return bScientific; }

    class column_iterator {
    public:
        using self_type = column_iterator;

        // conform to the iterator traits
        // https://en.cppreference.com/w/cpp/iterator/iterator_traits
        using difference_type = std::ptrdiff_t;
        using value_type = std::pair<REAL &, size_type const &>;
        using pointer = value_type *;
        using reference = value_type &;
        using iterator_category = std::bidirectional_iterator_tag;

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

        [[nodiscard]] value_type operator*() {
            return std::make_pair(std::ref(*_valIter),
                                  std::cref(*_colIndicesIter));
        }
        // [[nodiscard]] value_type operator->() {
        //     return std::make_pair(std::ref(*_valIter),
        //                           std::cref(*_colIndicesIter));
        // }

        [[nodiscard]] typename value_type::first_type value() {
            return std::ref(*_valIter);
        }

        [[nodiscard]] typename value_type::second_type index() {
            return std::cref(*_colIndicesIter);
        }

        [[nodiscard]] bool operator==(const self_type &other) {
            return (_valIter == other._valIter) and
                   (_colIndicesIter == other._colIndicesIter);
        }
        [[nodiscard]] bool operator!=(const self_type &other) {
            return not (*this == other);
        }

    private:
        VectorIterator _valIter;
        std::vector<size_type>::iterator _colIndicesIter;
    };

    class const_column_iterator {
    public:
        using self_type = const_column_iterator;

        // conform to the iterator traits
        // https://en.cppreference.com/w/cpp/iterator/iterator_traits
        using difference_type = std::ptrdiff_t;
        using value_type = std::pair<REAL const &, size_type const &>;
        using pointer = value_type *;
        using reference = value_type &;
        using iterator_category = std::bidirectional_iterator_tag;

        const_column_iterator(
            ConstVectorIterator valIter,
            std::vector<size_type>::const_iterator colIndicesIter)
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

        [[nodiscard]] value_type operator*() {
            return std::make_pair(std::ref(*_valIter),
                                  std::cref(*_colIndicesIter));
        }
        // TODO: This is wrong
        // [[nodiscard]] value_type operator->() {
        //     return std::make_pair(*_valIter, *_colIndicesIter);
        // }

        [[nodiscard]] typename value_type::first_type value() {
            return std::ref(*_valIter);
        }

        [[nodiscard]] typename value_type::second_type index() {
            return std::cref(*_colIndicesIter);
        }

        [[nodiscard]] bool operator==(const self_type &other) {
            return (_valIter == other._valIter) and
                   (_colIndicesIter == other._colIndicesIter);
        }
        [[nodiscard]] bool operator!=(const self_type &other) {
            return not (*this == other);
        }

    private:
        ConstVectorIterator _valIter;
        std::vector<size_type>::const_iterator _colIndicesIter;
    };

    class row_iterator {
    public:
        using self_type = row_iterator;

        // conform to the iterator traits
        // https://en.cppreference.com/w/cpp/iterator/iterator_traits
        using difference_type = std::ptrdiff_t;
        using value_type = self_type;
        using pointer = self_type *;
        using reference = self_type &;
        using iterator_category = std::random_access_iterator_tag;

        row_iterator(std::vector<size_type>::iterator rowPtrIter,
                     std::vector<size_type>::iterator colIndicesIter,
                     VectorIterator valIter)
            : _rowPtrIter(rowPtrIter), _colIndicesIter(colIndicesIter),
              _valIter(valIter) {}

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
            return *this;
        }

        // postfix
        self_type &operator++(int junk) {
            self_type cached = *this;
            _rowPtrIter++;
            return cached;
        }

        self_type &operator+=(difference_type offset) {
            _rowPtrIter += offset;
            return *this;
        }

        self_type &operator-=(difference_type offset) {
            _rowPtrIter -= offset;
            return *this;
        }

        // iter - n
        self_type operator-(difference_type offset) {
            self_type cache(*this);
            cache -= offset;
            return cache;
        }

        // iter + n
        self_type operator+(difference_type offset) {
            self_type cache(*this);
            cache += offset;
            return cache;
        }
        // n + iter
        friend self_type operator+(const difference_type &offset,
                                   const self_type &sec) {
            self_type cache(sec);
            cache += offset;
            return cache;
        }

        reference operator[](difference_type offset) {
            return *(*this + offset);
        }

        bool operator<(const self_type &other) {
            return other - (*this) > 0;  //
        }

        bool operator>(const self_type &other) {
            return other < (*this);  //
        }

        [[nodiscard]] self_type &operator*() { return *this; }
        // [[nodiscard]] self_type &operator->() { return *this; }

        [[nodiscard]] bool operator==(const self_type &rhs) {
            return _rowPtrIter == rhs._rowPtrIter;
        }
        [[nodiscard]] bool operator!=(const self_type &rhs) {
            return _rowPtrIter != rhs._rowPtrIter;
        }

    private:
        std::vector<size_type>::iterator _rowPtrIter;
        std::vector<size_type>::iterator _colIndicesIter;
        VectorIterator _valIter;
    };

    class const_row_iterator {
    public:
        using self_type = const_row_iterator;

        // conform to the iterator traits
        // https://en.cppreference.com/w/cpp/iterator/iterator_traits
        using difference_type = std::ptrdiff_t;
        using value_type = self_type;
        using pointer = self_type *;
        using reference = self_type &;
        using iterator_category = std::bidirectional_iterator_tag;

        const_row_iterator(
            std::vector<size_type>::const_iterator rowPtrIter,
            std::vector<size_type>::const_iterator colIndicesIter,
            ConstVectorIterator valIter)
            : _rowPtrIter(rowPtrIter), _colIndicesIter(colIndicesIter),
              _valIter(valIter) {}

        [[nodiscard]] const_column_iterator begin() const {
            return const_column_iterator((_valIter + *_rowPtrIter),
                                         (_colIndicesIter + *_rowPtrIter));
        }
        [[nodiscard]] const_column_iterator end() const {
            return const_column_iterator(
                (_valIter + *(_rowPtrIter + 1)),
                (_colIndicesIter + *(_rowPtrIter + 1)));
        }

        [[nodiscard]] const_column_iterator cbegin() const {
            return this->begin();
        }
        [[nodiscard]] const_column_iterator cend() const {
            return this->end();  //
        }

        // prefix
        self_type &operator++() {
            _rowPtrIter++;
            return *this;
        }

        // postfix
        self_type &operator++(int junk) {
            self_type cached = *this;
            _rowPtrIter++;
            return cached;
        }

        self_type &operator+=(difference_type offset) {
            _rowPtrIter += offset;
            return *this;
        }

        self_type &operator-=(difference_type offset) {
            _rowPtrIter -= offset;
            return *this;
        }

        // iter - n
        self_type operator-(difference_type offset) {
            self_type cache(*this);
            cache -= offset;
            return cache;
        }

        // iter + n
        self_type operator+(difference_type offset) {
            self_type cache(*this);
            cache += offset;
            return cache;
        }
        // n + iter
        friend self_type operator+(const difference_type &offset,
                                   const self_type &sec) {
            self_type cache(sec);
            cache += offset;
            return cache;
        }

        reference operator[](difference_type offset) {
            return *(*this + offset);
        }

        bool operator<(const self_type &other) {
            return other - (*this) > 0;  //
        }

        bool operator>(const self_type &other) {
            return other < (*this);  //
        }

        [[nodiscard]] self_type &operator*() { return *this; }
        // [[nodiscard]] self_type &operator->() { return this; }

        [[nodiscard]] bool operator==(const self_type &rhs) {
            return _rowPtrIter == rhs._rowPtrIter;
        }
        [[nodiscard]] bool operator!=(const self_type &rhs) {
            return _rowPtrIter != rhs._rowPtrIter;
        }

    private:
        std::vector<size_type>::const_iterator _rowPtrIter;
        std::vector<size_type>::const_iterator _colIndicesIter;
        ConstVectorIterator _valIter;
    };

    // regular (possibly modifying) Iterators
    [[nodiscard]] row_iterator begin() {
        return row_iterator(_rowPtr.begin(), _colIndices.begin(),
                            _data.begin());
    }
    [[nodiscard]] row_iterator end() {
        return row_iterator(_rowPtr.end() - 1, _colIndices.begin(),
                            _data.begin());
    }

    // const Iterators
    [[nodiscard]] const_row_iterator cbegin() const {
        return const_row_iterator(_rowPtr.cbegin(), _colIndices.cbegin(),
                                  _data.cbegin());
    }
    [[nodiscard]] const_row_iterator cend() const {
        return const_row_iterator(_rowPtr.cend() - 1, _colIndices.cbegin(),
                                  _data.cbegin());
    }
    [[nodiscard]] const_row_iterator begin() const { return this->cbegin(); }
    [[nodiscard]] const_row_iterator end() const { return this->cend(); }

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

    bool checkIfAccessIsInBounds(const size_type row_index,
                                 const size_type col_index) const {
        if (not (col_index < m_cols)) {
            HDNUM_ERROR("Out of bounds access: column too big! -> " +
                        std::to_string(col_index) + " is not < " +
                        std::to_string(m_cols));
            return false;
        } else if (not (row_index < m_rows)) {
            HDNUM_ERROR("Out of bounds access: row too big! -> " +
                        std::to_string(row_index) + " is not < " +
                        std::to_string(m_rows));
            return false;
        }
        return true;
    }

    // write access on matrix element A_ij using A.get(i,j)
    REAL &get(const size_type row_index, const size_type col_index) {
        checkIfAccessIsInBounds(row_index, col_index);
        // look for the entry
        using value_pair = typename const_column_iterator::value_type;
        auto row = row_iterator(_rowPtr.begin() + row_index,
                                _colIndices.begin(), _data.begin());
        auto result =
            std::find_if(row.begin(), row.end(), [col_index](value_pair el) {
                // only care for the index here
                // since the value is unknown
                // anyways
                return el.second == col_index;
            });
        // we found something within the right row
        if (result != row.end()) {
            return result.value();
        }
        throw std::out_of_range(
            "There is no non-zero element for these given indicies!");
    }

    //! read-access on matrix element A_ij using A(i,j)
    const REAL &operator()(const size_type row_index,
                           const size_type col_index) const {
        checkIfAccessIsInBounds(row_index, col_index);

        using value_pair = typename const_column_iterator::value_type;
        auto row = const_row_iterator(_rowPtr.begin() + row_index,
                                      _colIndices.begin(), _data.begin());
        auto result =
            std::find_if(row.begin(), row.end(), [col_index](value_pair el) {
                // only care for the index here since the value is unknown
                return el.second == col_index;
            });
        // we found something within the right row
        if (result != row.end()) {
            return result.value();
        }
        return _zero;
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
        SparseMatrix::builder builder(m_cols, m_rows);
        SparseMatrix::size_type curr_row = 0;
        for (auto &row : (*this)) {
            curr_row++;
            for (const std::pair<REAL const &, const size_type> pair : row) {
                builder.addEntry(pair.second, curr_row, pair.first);
            }
        }

        return builder.build();
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

    [[nodiscard]] std::string to_string() const noexcept {
        return "values=" + comma_fold(_data) + "\n" +        //
               "colInd=" + comma_fold(_colIndices) + "\n" +  //
               "rowPtr=" + comma_fold(_rowPtr) + "\n";       //
    }

    void print() const noexcept { std::cout << *this; }

    static SparseMatrix identity(const size_type dimN) {
        auto builder = typename SparseMatrix<REAL>::builder(dimN, dimN);
        for (typename SparseMatrix<REAL>::size_type i = 0; i < dimN; ++i) {
            builder.addEntry(i, i, REAL {1});
        }
        return builder.build();
    }

    SparseMatrix<REAL> matchingIdentity() const { return identity(m_cols); }

    class builder {
        size_type m_rows {};  // Number of Matrix rows, 0 by default
        size_type m_cols {};  // Number of Matrix columns, 0 by default
        std::vector<std::map<size_type, REAL>> _rows;

    public:
        builder(size_type new_m_rows, size_type new_m_cols)
            : m_rows {new_m_rows}, m_cols {new_m_cols}, _rows {m_rows} {}

        builder(const std::initializer_list<std::initializer_list<REAL>> &v)
            : m_rows {v.size()}, m_cols {v.begin()->size()}, _rows(m_rows) {
            size_type i = 0;
            for (auto &row : v) {
                size_type j = 0;
                for (const REAL &element : row) {
                    addEntry(i, j, element);
                    j++;
                }
                i++;
            }
        }

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
            return (m_rows == other.m_rows) and  //
                   (m_cols == other.m_cols) and  //
                   (_rows == other._rows);
        }

        [[nodiscard]] bool operator!=(
            const SparseMatrix::builder &other) const {
            return not (*this == other);
        }

        [[nodiscard]] size_type colsize() const noexcept { return m_cols; }
        [[nodiscard]] size_type rowsize() const noexcept { return m_rows; }

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

        [[nodiscard]] std::string to_string() const {
            std::string output;
            for (std::size_t i = 0; i < _rows.size(); i++) {
                for (const auto &indexpair : _rows[i]) {
                    output += std::to_string(i) + ", " +
                              std::to_string(indexpair.first) + " => " +
                              std::to_string(indexpair.second) + "\n";
                }
            }
            return output;
        }

        [[nodiscard]] SparseMatrix build() {
            auto result = SparseMatrix<REAL>(m_rows, m_cols);

            for (std::size_t i = 0; i < _rows.size(); i++) {
                result._rowPtr[i + 1] = result._rowPtr[i];
                for (const auto &indexpair : _rows[i]) {
                    result._colIndices.push_back(indexpair.first);
                    result._data.push_back(indexpair.second);
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
const REAL SparseMatrix<REAL>::_zero {};

template <typename REAL>
std::ostream &operator<<(std::ostream &s, const SparseMatrix<REAL> &A) {
    s << std::endl;
    s << " " << std::setw(A.iwidth()) << " "
      << "  ";
    for (typename SparseMatrix<REAL>::size_type j = 0; j < A.colsize(); ++j)
        s << std::setw(A.width()) << j << " ";
    s << std::endl;

    for (typename SparseMatrix<REAL>::size_type i = 0; i < A.rowsize(); ++i) {
        s << " " << std::setw(A.iwidth()) << i << "  ";
        for (typename SparseMatrix<REAL>::size_type j = 0; j < A.colsize();
             ++j) {
            if (A.scientific()) {
                s << std::setw(A.width()) << std::scientific << std::showpoint
                  << std::setprecision(A.precision()) << A(i, j) << " ";
            } else {
                s << std::setw(A.width()) << std::fixed << std::showpoint
                  << std::setprecision(A.precision()) << A(i, j) << " ";
            }
        }
        s << std::endl;
    }
    return s;
}

//! make a zero matrix
template <typename REAL>
inline void zero(SparseMatrix<REAL> &A) {
    A = SparseMatrix<REAL>();
}

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
template <class REAL>
inline void identity(SparseMatrix<REAL> &A) {
    if (A.rowsize() != A.colsize()) {
        HDNUM_ERROR("Will not overwrite A since Dimensions are not equal!");
    }
    A = SparseMatrix<REAL>::identity(A.colsize());
}

template <typename REAL>
inline void readMatrixFromFile(const std::string &filename,
                               SparseMatrix<REAL> &A) {
    // Format taken from here:
    // https://math.nist.gov/MatrixMarket/formats.html#coord

    using size_type = typename SparseMatrix<REAL>::size_type;
    std::string buffer;
    std::ifstream fin(filename);
    size_type i = 0;
    size_type j = 0;
    size_type non_zeros = 0;

    if (fin.is_open()) {
        // ignore all comments from the file (starting with %)
        while (fin.peek() == '%') fin.ignore(2048, '\n');

        std::getline(fin, buffer);
        std::istringstream first_line(buffer);
        first_line >> i >> j >> non_zeros;

        auto builder = typename SparseMatrix<REAL>::builder(i, j);

        while (std::getline(fin, buffer)) {
            std::istringstream iss(buffer);

            REAL value {};
            iss >> i >> j >> value;

            builder.addEntry(i - 1, j - 1,
                             value);  // i-1, j-1, because matrix market does
                                      // not use zero based indexing
        }
        A = builder.build();
        fin.close();
    } else {
        HDNUM_ERROR("Could not open file!");
    }
}

}  // namespace hdnum

#endif  // SPARSEMATRIX_HH
