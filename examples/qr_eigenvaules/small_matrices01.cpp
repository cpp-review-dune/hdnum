
#include <iostream>    
#include "../../../hdnum/src/qr_eigenvalues.hh" //TODO change to cmake


using namespace hdnum;

using Number = double;
using Mat = DenseMatrix<Number>;


int main(){
    Mat A= { {1,2,3 }, {4, 5, 6}, {7, 8, 9}};
    Mat B= { {1,2,3, 3 }, {4, 5, 6, 6}, {7, 8, 9, 9}, {1,1,14,1}};
    Mat C= { {1.0,2,3, 3, 2 }, {3, 4, 5, 6, 6}, {7, 8,6.7, 9, 9}, {1, 5,12,14,1}, { 34, 1, 6, 9, 5}};

    std::cout <<C << std::endl;


    QR_Info<Number> qr_info = eigenvalues_qr_algorithm_givens(C);
    std::vector<Number> real = qr_info.get_real();
    std::vector<Number> imag = qr_info.get_imag();

    for (int i = 0; i< real.size(); i++){
        std::cout << real[i]  << " + i*" << imag[i] << std::endl;
    }
    return 0;
}