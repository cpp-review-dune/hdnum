#include <vector>
#include "../../hdnum.hh"
#include "laplace_matrix.hh"

using namespace hdnum;

int main ()
{
  typedef double Number;
  size_t n = 3; // DoFs per dimension

  auto M = hdnum::laplaceMatrixDimSparse<Number>(n,6);
  std::cout << M;

  return 0;
}
