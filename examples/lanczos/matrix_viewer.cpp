#include <iostream>
#include "../hdnum.hh"
#include <string>

/*
 * symmetric matrices:
   bcsstm01
   bcsstm22
   Alemdar

 * 
 */

using namespace hdnum;

int main(){
   
    std::string path = "../test/matrix_market_files/ex_lanczos_nos4.mtx";

    SparseMatrix<double> sampleMatrix {};
    readMatrixFromFile(path, sampleMatrix);
    bool is_symmetric = true;

    for(long i; i < sampleMatrix.colsize(); i++) {
        for(long j; j < sampleMatrix.colsize(); j++){
            if(sampleMatrix(i,j) != sampleMatrix(j,i)) {
                is_symmetric = false;
            }
        }
    }

    if (is_symmetric) {
        std::cout << "Matrix is symmetric!" << std::endl;
    }
    else {
        std::cout << "Matrix is not symmetric!" << std::endl;
    }

    sampleMatrix.scientific(false);
    sampleMatrix.width(1);
    sampleMatrix.precision(1);
    //std::cout << "SampleMatrix: \n" << sampleMatrix << std::endl;

    return 0;
}