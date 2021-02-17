#include <algorithm>
#include <utility>
#include <vector>

#include "../hdnum.hh"
#include "gtest/gtest.h"

namespace {
// In this example, we test the operators for the SparseMatrix class.

using hdnum::SparseMatrix;

class TestReadSparseMatrixFromFile : public ::testing::Test {
public:
    const SparseMatrix<double> A;
    const std::string filename = "../example_matrix_market.mtx";

    const SparseMatrix<double> A_empty;
    const std::string filename_empty = "../example_empty_matrix_market.mtx";

    using size_type = typename SparseMatrix<double>::size_type;
};

TEST_F(TestReadSparseMatrixFromFile, ReadInMatrixFromEmptyFile) {
    SparseMatrix<double> B_empty {};
    hdnum::readMatrixFromFile(this->filename_empty, B_empty);

    // check sparse matrix is empty
    ASSERT_EQ(B_empty.rowsize(), 0);
    ASSERT_EQ(B_empty.colsize(), 0);
}

TEST_F(TestReadSparseMatrixFromFile, ReadInMatrixFromValidFile) {
    SparseMatrix<double> B {};
    hdnum::readMatrixFromFile(this->filename, B);

    ASSERT_EQ(B.rowsize(), 5);
    ASSERT_EQ(B.colsize(), 5);

    std::vector<std::pair<int, int>> nonzeros = {
        {0, 0}, {4, 0}, {1, 1}, {2, 2}, {3, 3}};

    // check nonzero elements
    ASSERT_EQ(B(0, 0), 1474.779);
    ASSERT_EQ(B(4, 0), -9.017133);
    ASSERT_EQ(B(1, 1), 9.136654);
    ASSERT_EQ(B(2, 2), 69.61468);
    ASSERT_EQ(B(3, 3), 68.60106);

    // check zero elements
    for (int r = 0; r < B.rowsize(); r++) {
        for (int c = 0; c < B.colsize(); c++) {
            // nonzero
            if (std::find(nonzeros.begin(), nonzeros.end(),
                          std::make_pair(r, c)) != nonzeros.end()) {
                continue;
            }
            // zero
            EXPECT_THROW(
                {
                    try {
                        B(r, c);
                    } catch (const hdnum::ErrorException& e) {
                        throw;
                    }
                },
                std::out_of_range);
        }
    }
}

}  // namespace
