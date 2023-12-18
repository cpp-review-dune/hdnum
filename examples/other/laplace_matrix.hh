// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/*
 * File:   laplace_matrix.hh
 * Author: Michal Tóth
 *
 * Created on December 18, 2023
 */

#ifndef LAPLACEMATRIX_HH
#define LAPLACEMATRIX_HH

namespace hdnum {

  /*! \brief Create a Laplace matrix of size n^2
   *
   * Laplace matrix for 2D Poisson problem
   * with n degrees of freedom in each direction.
   *
   * Matrix size n^2
   * Diagonal entries of value 4
   * 4 offdiagonals with entries -1
   */
  template<typename Number>
  DenseMatrix<Number> laplaceMatrix (size_t n)
  {
    DenseMatrix<Number> M(n*n,n*n,0);
    // iterate over blocks
    for (size_t i=0; i<n; ++i)
    {
      // iterate inside the block
      for (size_t j=0; j<n; ++j)
      {
        M(i*n+j, i*n+j) = 4;
        if (j > 0  )  M(i*n+j, i*n+j-1) = -1; // sub-diagonal
        if (j < n-1)  M(i*n+j, i*n+j+1) = -1; // over-diagonal
        if (i > 0  )  M(i*n+j, i*n+j-n) = -1; // block sub-diagonal
        if (i < n-1)  M(i*n+j, i*n+j+n) = -1; // block over-diagonal
      }
    }
    return M;
  }

  template<typename Number>
  SparseMatrix<Number> laplaceMatrixSparse (size_t n)
  {
    typename SparseMatrix<Number>::builder builder(n*n, n*n);
    // iterate over blocks
    for (size_t i=0; i<n; ++i)
    {
      // iterate inside the block
      for (size_t j=0; j<n; ++j)
      {
        builder.addEntry(i*n+j, i*n+j, 4);
        if (j > 0  )  builder.addEntry(i*n+j, i*n+j-1, -1); // sub-diagonal
        if (j < n-1)  builder.addEntry(i*n+j, i*n+j+1, -1); // over-diagonal
        if (i > 0  )  builder.addEntry(i*n+j, i*n+j-n, -1); // block sub-diagonal
        if (i < n-1)  builder.addEntry(i*n+j, i*n+j+n, -1); // block over-diagonal
      }
    }
    return builder.build();
  }

  /*! \brief Create a Laplace matrix of arbitrary dimension of size n^dim
   *
   * Laplace matrix for Poisson problem on a dim-dimensional cube
   * with n degrees of freedom in each direction.
   *
   * Matrix size n^dim
   * Diagonal entries of value 2*dim
   * 2*dim offdiagonals with entries -1
   */
  template<typename Number>
  DenseMatrix<Number> laplaceMatrixDim (size_t n, size_t dim=2)
  {
    size_t size = pow(n,dim);
    DenseMatrix<Number> M(size,size,0);

    std::vector<size_t> powers; // offdiagonals are n^i away from diagonal, i=0..dim-1
    for (size_t i=1; i<size; i*=n)
      powers.push_back(i);

    for (size_t i=0; i<size; ++i)
    {
      M(i,i) = 2*dim; // diagonal

      for (const auto& v : powers)
      {
        // offdiagonals have n^dim-n^(dim-1) entries in the block structure,
        // omitted entries correspond to DoFs on the cube's face
        if (i   >= v   && (i/v)%n != 0  )  M(i, i-v) = -1;
        if (i+v < size && (i/v)%n != n-1)  M(i, i+v) = -1;
      }
    }

    return M;
  }

  template<typename Number>
  SparseMatrix<Number> laplaceMatrixDimSparse (size_t n, size_t dim=2)
  {
    size_t size = pow(n,dim);
    typename SparseMatrix<Number>::builder builder(size, size);

    std::vector<size_t> powers; // offdiagonals are n^i away from diagonal, i=0..dim-1
    for (size_t i=1; i<size; i*=n)
      powers.push_back(i);

    for (size_t i=0; i<size; ++i)
    {
      builder.addEntry(i,i, 2*dim); // diagonal

      for (const auto& v : powers)
      {
        // offdiagonals have n^dim-n^(dim-1) entries in the block structure,
        // omitted entries correspond to DoFs on the cube's face
        if (i   >= v   && (i/v)%n != 0  )  builder.addEntry(i, i-v, -1);
        if (i+v < size && (i/v)%n != n-1)  builder.addEntry(i, i+v, -1);
      }
    }

    return builder.build();
  }

} // end namespace hdnum

#endif // LAPLACEMATRIX_HH