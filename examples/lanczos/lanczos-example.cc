#include <iostream>
#include "../../hdnum.hh"
#include <string>

using namespace hdnum;

void eig(SparseMatrix<double> &Tk) {
    std::vector<double> real;
    std::vector<double> imag;
    DenseMatrix<double> Tk_dense = Tk;
    eigenvalues_qr_algorithm_givens(Tk_dense, real, imag);

    std::cout << "Eigenvalues: " << std::endl;
    for (int i = 0; i< real.size(); i++){
        std::cout << "[" << i << "]: "<< real[i]  << " + i*" << imag[i] << std::endl;
    } 
}

void matrix_output(SparseMatrix<double> &Tk) {
    Tk.scientific(false);
    Tk.width(8);
    Tk.precision(3);
    if (Tk.colsize() < 40) {
        std::cout << Tk << std::endl;
    }
}

void sample_run(SparseMatrix<double> sampleMatrix) {
    SparseMatrix<double> Tk_cro = lanczos_cro(sampleMatrix);

    sampleMatrix.scientific(false);
    sampleMatrix.width(8);
    sampleMatrix.precision(3);
    if (sampleMatrix.colsize() < 40){
        std::cout << "SampleMatrix: \n" << sampleMatrix << std::endl;
    }
    else {
        std::cout << "SampleMatrix with " << sampleMatrix.colsize() << "x" << sampleMatrix.rowsize() << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;

    // //Optimized Lanczos method
    SparseMatrix<double> Tk_opt = lanczos_optimized(sampleMatrix);
    std::cout << "Optimized Lanczos: \n" << std::endl;
    //matrix_output(Tk_opt);
    eig(Tk_opt);
    std::cout << "--------------------------------------------" << std::endl;

    // //Basic Lanczos method
    SparseMatrix<double> Tk_basic = lanczos_basic(sampleMatrix);
    std::cout << "Basic Lanczos: \n" << std::endl;
    //matrix_output(Tk_basic);
    eig(Tk_basic);
    std::cout << "--------------------------------------------" << std::endl;

    // Lanczos complete reorthogonalization method
    //SparseMatrix<double> Tk_cro = lanczos_cro(sampleMatrix);
    std::cout << "CRO Lanczos: \n" << std::endl;
    //matrix_output(Tk_cro);
    eig(Tk_cro);
    std::cout << "--------------------------------------------" << std::endl;
}

void sampleMatrix1(){
    //sample from: https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample1_simple.cpp
    std::initializer_list<std::initializer_list<double>> sample = {{2, 1, 1},
                                                                   {1, 2, 1},
                                                                   {1, 1, 2}};
    // Its eigenvalues are {4, 1, 1}
    SparseMatrix<double> sampleMatrix = SparseMatrix<double>::builder(sample).build();
    sample_run(sampleMatrix);
}

void sampleMatrix2(){
    std::initializer_list<std::initializer_list<double>> sample = {{0, 1, 1},
                                                                   {1, 0, -1},
                                                                   {1, -1, 0}};
    // Its eigenvalues are {1, 1, -2}
    SparseMatrix<double> sampleMatrix = SparseMatrix<double>::builder(sample).build();
    sample_run(sampleMatrix);
}

void sampleMatrixFromCollection(std::string &path){
    SparseMatrix<double> sampleMatrix {};
    readMatrixFromFile(path, sampleMatrix);
    sample_run(sampleMatrix);
}

int main(){
    std::string path;
    //sampleMatrix1 (3x3):
    sampleMatrix1();
    // Its eigenvalues are {4, 1, 1}

    //sampleMatrix2 (3x3):
    sampleMatrix2();
    // Its eigenvalues are {1, 1, -2}

    //Suite Sparse Matrix Collection
    //sampleMatrix1 bcsstm01 (48x48):
    path = "matrix_files/ex_lanczos_bcsstm01.mtx";
    // Its eigenvalues are {200, 200, 200, 200, 200, 200, ..., 100, ..., 0}
    //sampleMatrixFromCollection(path);

    //sampleMatrix2 bcsstm22 (138x138):
    path = "matrix_files/ex_lanczos_bcsstm22.mtx";
    // Its eigenvalues are {0.00973766582616700, 0.00973766582616699, 0.00973766582616601, 0.00973766582616600, 0.00973766582616599, 0.00973766582616598, 0.00835385206984501}
    //sampleMatrixFromCollection(path);

    //sampleMatrix3 Alemdar (6245x6245):
    path = "matrix_files/ex_lanczos_Alemdar.mtx";
    // Its eigenvalues are {69.5187762679672, 69.2682886007168, 68.9967100511161, 68.9689724618314, 68.7323042832841, 68.6419265794814, 68.5722507550780}
    //sampleMatrixFromCollection(path);

    //sampleMatrix12 aug3d (24300x24300):
    path = "matrix_files/ex_lanczos_aug3d.mtx";
    // Its eigenvalues are {3.45226995184565, -3.45226995184564, 3.44050573273934, 3.44050573273934, 3.44050573273933, -3.44050573273933, -3.44050573273932}
    //sampleMatrixFromCollection(path);

    return 0;
}