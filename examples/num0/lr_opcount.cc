#include <iostream>
#include "hdnum.hh"

int main ()
{
  // Setup matrix for LU decomposition
  const int n=7;
  typedef hdnum::oc::OpCounter<double> number;
  hdnum::DenseMatrix<number> A(n,n);
  hdnum::spd(A);
  hdnum::Vector<std::size_t> perm(n);

  // Make LU decomposition and store number of floating point operations
  std::ostringstream os;
  A[0][0].reset();
  hdnum::lr(A, perm);
  A[0][0].reportOperations(os);
  int operations = A[0][0].totalOperationCount();

  // Print number of operations
  std::cout << "Total number of operations: " << operations << std::endl;
  std::cout << "In detail:" << std::endl << os.str() << std::endl;
}
