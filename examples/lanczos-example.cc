#include <iostream>
#include "../hdnum.hh"
#include <string>

using namespace hdnum;

void sampleRun(SparseMatrix<double> sampleMatrix) {
    sampleMatrix.scientific(false);
    sampleMatrix.width(8);
    sampleMatrix.precision(4);
    std::cout << "SampleMatrix: \n" << sampleMatrix << std::endl;

    SparseMatrix<double> Tk = lanczos(sampleMatrix);
    Tk.scientific(false);
    Tk.width(8);
    Tk.precision(4);
    std::cout << "Optimized Lanczos: \n" << Tk << std::endl;

//     SparseMatrix<double> Tk_basic = lanczos_basic(sampleMatrix);
//     Tk_basic.scientific(false);
//     Tk_basic.width(8);
//     Tk_basic.precision(4);
//     std::cout << "Basic Lanczos: \n" << Tk_basic << std::endl;
//     std::cout << "--------------------------------------------" << std::endl;
}

void sampleMatrix1(){
    //sample from: https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample1_simple.cpp
    std::initializer_list<std::initializer_list<double>> sample = {{2, 1, 1},
                                                                   {1, 2, 1},
                                                                   {1, 1, 2}};
    // Its eigenvalues are {4, 1, 1}
    SparseMatrix<double> sampleMatrix = SparseMatrix<double>::builder(sample).build();
    sampleRun(sampleMatrix);
}

void sampleMatrix2(){
    std::initializer_list<std::initializer_list<double>> sample = {{0, 1, 1},
                                                                   {1, 0, -1},
                                                                   {1, -1, 0}};
    // Its eigenvalues are {1, 1, -2}
    SparseMatrix<double> sampleMatrix = SparseMatrix<double>::builder(sample).build();
    sampleRun(sampleMatrix);
}

void sampleMatrixFromCollection(std::string &path){
    SparseMatrix<double> sampleMatrix {};
    readMatrixFromFile(path, sampleMatrix);
    sampleRun(sampleMatrix);
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
    //sampleMatrix3 bcsstm01 (48x48):
    path = "../test/matrix_market_files/ex_lanczos_bcsstm01.mtx";
    // Its eigenvalues are {100, 200, 0}
    //sampleMatrixFromCollection(path);

    //sampleMatrix11 mycielskian3 (5x5):
    path = "../test/matrix_market_files/ex_lanczos_mycielskian3.mtx";
    // Its eigenvalues are {}
    //sampleMatrixFromCollection(path);

    //sampleMatrix12 Trefethen_20b (19x19):
    path = "../test/matrix_market_files/ex_lanczos_Trefethen_20b.mtx";
    // Its eigenvalues are {}
    sampleMatrixFromCollection(path);

    //sampleMstrix6 bcsstk01 (48x48):
    path = "../test/matrix_market_files/ex_lanczos_bcsstk01.mtx";
    //sampleMatrixFromCollection(path);

    //sampleMatrix5 bcspwr02 (49x49):
    path = "../test/matrix_market_files/ex_lanczos_bcspwr02.mtx";
    //sampleMatrixFromCollection(path);

    //sampleMatrix4 bcsstm22 (138x138):
    path = "../test/matrix_market_files/ex_lanczos_bcsstm22.mtx";
    // Its eigenvalues are {}
    //sampleMatrixFromCollection(path);

    //sampleMatrix7 bcsstk07 (420x420):
    path = "../test/matrix_market_files/ex_lanczos_bcsstk07.mtx";
    //sampleMatrixFromCollection(path);

    //sampleMatrix10 G11 (800x800):
    path = "../test/matrix_market_files/ex_lanczos_G11.mtx";
    // Its eigenvalues are {}
    //sampleMatrixFromCollection(path);

    //sampleMatrix8 bcsstk13 (2003x2003):
    path = "../test/matrix_market_files/ex_lanczos_bcsstk13.mtx";
    //sampleMatrixFromCollection(path);

    //sampleMatrix9 ex10 (2410x2410):
    path = "../test/matrix_market_files/ex_lanczos_ex10.mtx";
    // Its eigenvalues are {}
    //sampleMatrixFromCollection(path);
   
    return 0;
}