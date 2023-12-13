#include <iostream>
#include "../../hdnum.hh"
#include <string>

using namespace hdnum;

int main(){
    // std::initializer_list<std::initializer_list<double>> sample = {{1, 2, 4},
    //                                                                {0, 0, 5},
    //                                                                {0, 3, 6}};
    // SparseMatrix<double> sampleMatrix = SparseMatrix<double>::builder(sample).build();

    std::initializer_list<std::initializer_list<double>> sample = {{1, -4},
                                                                   {2, 3},
                                                                   {2, 2}};
    SparseMatrix<double> sampleMatrix = SparseMatrix<double>::builder(sample).build();

    sampleMatrix.scientific(false);
    sampleMatrix.width(8);
    sampleMatrix.precision(3);
    std::cout << "SampleMatrix: \n" << sampleMatrix << std::endl;
    
    Vector<double> x(3);
    x[0] = 1; x[1] = 2; x[2] = 2;
    std::cout << "x: " << x << std::endl;

    Vector<double> estimated_v(3);
    estimated_v[0] = 2; estimated_v[1] = 0; estimated_v[2] = 0;
    std::cout << "Estimated v: " << estimated_v << std::endl;

    std::pair<Vector<double>, double> hh = householder(x);

    std::cout << "Computed v: " << hh.first << std::endl;

    DenseMatrix<double> H (3,3, 0);
    DenseMatrix<double> I(3,3);
    identity(I);

    DenseMatrix<double> v_matrix(3,1);
    for (int i = 0; i < 3; i++) {
        v_matrix[i][1] = hh.first[i];
    }
    DenseMatrix<double> vT(1,3);
    vT = v_matrix.transpose();

    DenseMatrix<double> vTv(3,3);
    vTv.mm(v_matrix, vT);
    vTv *= hh.second;

    H = I;
    H -= vTv;

    std::cout << "H: " << H << std::endl;

    return 0;
}