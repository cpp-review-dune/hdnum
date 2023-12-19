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
   * Laplace matrix for 2D Poisson problem on a square with Dirichlet
   * boundary conditions and with n degrees of freedom in each direction
   * (n+2 points when counting boundary points which do not enter the matrix).
   *
   * Matrix size n^2
   * Diagonal entries of value 4
   * 4 offdiagonals with entries -1
   * Can rescale to unit square by adding the factor (n+1)^2
   */
  template<typename Number>
  DenseMatrix<Number> laplaceMatrix (size_t n, bool rescaled=true)
  {
    DenseMatrix<Number> M(n*n,n*n,0);
    Number inv_h2 = 1.;
    if (!rescaled) inv_h2 = (n+1)*(n+1); // inverse cell width squared
    // iterate over blocks
    for (size_t i=0; i<n; ++i)
    {
      // iterate inside the block
      for (size_t j=0; j<n; ++j)
      {
        M(i*n+j, i*n+j) = 4*inv_h2;
        if (j > 0  )  M(i*n+j, i*n+j-1) = -1*inv_h2; // sub-diagonal
        if (j < n-1)  M(i*n+j, i*n+j+1) = -1*inv_h2; // over-diagonal
        if (i > 0  )  M(i*n+j, i*n+j-n) = -1*inv_h2; // block sub-diagonal
        if (i < n-1)  M(i*n+j, i*n+j+n) = -1*inv_h2; // block over-diagonal
      }
    }
    return M;
  }

  template<typename Number>
  SparseMatrix<Number> laplaceMatrixSparse (size_t n, bool rescaled=true)
  {
    typename SparseMatrix<Number>::builder builder(n*n, n*n);
    Number inv_h2 = 1.;
    if (!rescaled) inv_h2 = (n+1)*(n+1); // inverse cell width squared
    // iterate over blocks
    for (size_t i=0; i<n; ++i)
    {
      // iterate inside the block
      for (size_t j=0; j<n; ++j)
      {
        builder.addEntry(i*n+j, i*n+j, 4*inv_h2);
        if (j > 0  )  builder.addEntry(i*n+j, i*n+j-1, -1*inv_h2); // sub-diagonal
        if (j < n-1)  builder.addEntry(i*n+j, i*n+j+1, -1*inv_h2); // over-diagonal
        if (i > 0  )  builder.addEntry(i*n+j, i*n+j-n, -1*inv_h2); // block sub-diagonal
        if (i < n-1)  builder.addEntry(i*n+j, i*n+j+n, -1*inv_h2); // block over-diagonal
      }
    }
    return builder.build();
  }

  /*! \brief Create a Laplace matrix of arbitrary dimension of size n^dim
   *
   * Laplace matrix for Poisson problem on a dim-dimensional unit cube with Dirichlet
   * boundary conditions and with n degrees of freedom in each direction
   * (n+2 points when counting boundary points which do not enter the matrix).
   *
   * Matrix size n^dim
   * Diagonal entries of value 2*dim
   * 2*dim offdiagonals with entries -1
   * Can rescale to unit square by adding the factor (n+1)^2
   */
  template<typename Number>
  DenseMatrix<Number> laplaceMatrixDim (size_t n, size_t dim=2, bool rescaled = true)
  {
    size_t size = pow(n,dim);
    DenseMatrix<Number> M(size,size,0);
    Number inv_h2 = 1.;
    if (!rescaled) inv_h2 = (n+1)*(n+1); // inverse cell width squared

    std::vector<size_t> powers; // offdiagonals are n^i away from diagonal, i=0..dim-1
    for (size_t i=1; i<size; i*=n)
      powers.push_back(i);

    for (size_t i=0; i<size; ++i)
    {
      M(i,i) = 2*dim*inv_h2; // diagonal

      for (const auto& v : powers)
      {
        // offdiagonals have n^dim-n^(dim-1) entries in the block structure,
        // omitted entries correspond to DoFs on the cube's face
        if (i   >= v   && (i/v)%n != 0  )  M(i, i-v) = -1*inv_h2;
        if (i+v < size && (i/v)%n != n-1)  M(i, i+v) = -1*inv_h2;
      }
    }

    return M;
  }

  template<typename Number>
  SparseMatrix<Number> laplaceMatrixDimSparse (size_t n, size_t dim=2, bool rescaled=true)
  {
    size_t size = pow(n,dim);
    typename SparseMatrix<Number>::builder builder(size, size);
    Number inv_h2 = 1.;
    if (!rescaled) inv_h2 = (n+1)*(n+1); // inverse cell width squared

    std::vector<size_t> powers; // offdiagonals are n^i away from diagonal, i=0..dim-1
    for (size_t i=1; i<size; i*=n)
      powers.push_back(i);

    for (size_t i=0; i<size; ++i)
    {
      builder.addEntry(i,i, 2*dim*inv_h2); // diagonal

      for (const auto& v : powers)
      {
        // offdiagonals have n^dim-n^(dim-1) entries in the block structure,
        // omitted entries correspond to DoFs on the cube's face
        if (i   >= v   && (i/v)%n != 0  )  builder.addEntry(i, i-v, -1*inv_h2);
        if (i+v < size && (i/v)%n != n-1)  builder.addEntry(i, i+v, -1*inv_h2);
      }
    }

    return builder.build();
  }

  /*! \brief Evaluate eigenvalues of Laplace matrix
   *
   * Eigenvalues of matrix corresponding to Poisson problem on a cube with Dirichlet
   * boundary conditions and with n degrees of freedom in each direction
   * (n+2 points when counting boundary points which do not enter the matrix).
   *
   * n^dim eigenvalues sorted from lowest to highest.
   * Eigenmvalues correspond to rescaled matrix that has -1 on offdiagonals,
   * can be rescaled back to the unit square problem's eigenvalues.
   */
  template<typename Number>
  Vector<Number> eigenvaluesLaplace(size_t n, size_t dim=2, bool rescaled=true)
  {
    assert(n>0 && dim>0);

    // evaluate 1d eigenvalues, reused for dim>1
    Vector<Number> ev1d(n);
    Number h = 1 / (n + 1.0);
    for (std::size_t i = 0; i < n; ++i)
    {
      ev1d[i] = 2 - 2 * std::cos(h * (i+1) * M_PI);
      if (!rescaled) ev1d[i] *= (n+1)*(n+1);
    }

    if (dim==1)
    {
      std::sort(ev1d.begin(),ev1d.end());
      return ev1d;
    }

    size_t size = std::pow(n,dim); // number of eigenvalues
    Vector<Number> ev(size);
    // Problem is separable, in higher dimensions we just add eigenvalues
    // of 1D problem dim-times. All combinations yield all eigenvalues.
    // For example in 3D: ev[c*n*n + b*n +a] = ev1d[c] + ev1d[b] + ev1d[a];
    for (size_t i=0; i<size; ++i)
    {
      ev[i] = 0.;
      size_t k=i;
      for (size_t j=0; j<dim; ++j)
      {
        ev[i] += ev1d[k % n];
        k /= n;
      }
    }

    std::sort(ev.begin(), ev.end());
    return ev;
  }

} // end namespace hdnum

#endif // LAPLACEMATRIX_HH