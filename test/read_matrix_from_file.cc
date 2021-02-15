#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the operators for the SparseMatrix class.

using hdnum::SparseMatrix;

class TestReadSparseMatrixFromFile : public ::testing::Test {
public:
    const SparseMatrix<double> A;
    const std::string filename = "../example_matrix_market.mtx";

    using size_type = typename SparseMatrix<double>::size_type;
};

TEST_F(TestReadSparseMatrixFromFile, ReadInMatrixFromFile) {
    SparseMatrix<double> B {};
    hdnum::readMatrixFromFile(this->filename, B);
}
}  // namespace
