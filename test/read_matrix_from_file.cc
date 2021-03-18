#include <algorithm>
#include <utility>
#include <vector>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the read Matrix function from densematrix and
// sparsematrix

using hdnum::DenseMatrix;
using hdnum::SparseMatrix;

class TestReadMatrixFromFile : public ::testing::Test {
public:
    const SparseMatrix<double> A_sparse;
    const std::string filename =
        "../matrix_market_files/example_matrix_market.mtx";

    const DenseMatrix<double> A_dense;

    const SparseMatrix<double> A_sparse_empty;
    const std::string filename_empty =
        "../matrix_market_files/example_empty_matrix_market.mtx";

    const DenseMatrix<double> A_dense_empty;

    using size_type = typename SparseMatrix<double>::size_type;
};

TEST_F(TestReadMatrixFromFile, ReadInMatrixFromEmptyFile) {
    SparseMatrix<double> B_sparse_empty {};
    hdnum::readMatrixFromFile(this->filename_empty, B_sparse_empty);

    // check sparse matrix is empty
    ASSERT_EQ(B_sparse_empty.rowsize(), 0);
    ASSERT_EQ(B_sparse_empty.colsize(), 0);

    DenseMatrix<double> B_dense_empty {};
    hdnum::readMatrixFromFileMatrixMarket(this->filename_empty, B_dense_empty);

    // check sparse matrix is empty
    ASSERT_EQ(B_dense_empty.rowsize(), 0);
    ASSERT_EQ(B_dense_empty.colsize(), 0);
}

TEST_F(TestReadMatrixFromFile, ReadInMatrixFromValidFile) {
    SparseMatrix<double> B_sparse {};
    hdnum::readMatrixFromFile(this->filename, B_sparse);

    ASSERT_EQ(B_sparse.rowsize(), 5);
    ASSERT_EQ(B_sparse.colsize(), 5);

    std::vector<std::pair<int, int>> nonzeros = {
        {0, 0}, {4, 0}, {1, 1}, {2, 2}, {3, 3}};

    // check nonzero elements
    ASSERT_EQ(B_sparse(0, 0), 1474.779);
    ASSERT_EQ(B_sparse(4, 0), -9.017133);
    ASSERT_EQ(B_sparse(1, 1), 9.136654);
    ASSERT_EQ(B_sparse(2, 2), 69.61468);
    ASSERT_EQ(B_sparse(3, 3), 68.60106);

    // check zero elements
    for (int r = 0; r < B_sparse.rowsize(); r++) {
        for (int c = 0; c < B_sparse.colsize(); c++) {
            // nonzero
            if (std::find(nonzeros.begin(), nonzeros.end(),
                          std::make_pair(r, c)) != nonzeros.end()) {
                continue;
            }
            // zero
            ASSERT_EQ(B_sparse(r, c), 0);
        }
    }

    DenseMatrix<double> B_dense {};
    hdnum::readMatrixFromFileMatrixMarket(this->filename, B_dense);

    ASSERT_EQ(B_dense.rowsize(), 5);
    ASSERT_EQ(B_dense.colsize(), 5);

    // check nonzero elements
    ASSERT_EQ(B_dense(0, 0), 1474.779);
    ASSERT_EQ(B_dense(4, 0), -9.017133);
    ASSERT_EQ(B_dense(1, 1), 9.136654);
    ASSERT_EQ(B_dense(2, 2), 69.61468);
    ASSERT_EQ(B_dense(3, 3), 68.60106);

    // check zero elements
    for (int r = 0; r < B_dense.rowsize(); r++) {
        for (int c = 0; c < B_dense.colsize(); c++) {
            // nonzero
            if (std::find(nonzeros.begin(), nonzeros.end(),
                          std::make_pair(r, c)) != nonzeros.end()) {
                continue;
            }
            // zero
            ASSERT_EQ(B_dense(r, c), 0);
        }
    }
}

}  // namespace
