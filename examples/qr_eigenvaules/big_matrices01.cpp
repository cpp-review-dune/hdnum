
#include <iostream>    
#include "../../../hdnum/src/qr_eigenvalues.hh" //TODO change to cmake
#include "../../../hdnum/src/timer.hh" //TODO change to cmake



using namespace hdnum;


using Number = double;
using Mat = DenseMatrix<Number>;

Mat makeBigMatrix(int n){
    std::srand(1); //seed can be changed
    Mat A= Mat(n, n, 1);
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            A[i][j]=std::rand()/((RAND_MAX + 1u)/6);
        }
    }
    return A;
}

int main(){

    Mat A= makeBigMatrix(10); //TODO Error at 62
    A.scientific(false);


    std::vector<Number> real;
    std::vector<Number> imag;
    Timer t;
    t.reset();
    eigenvalues_qr_algorithm_givens(A, real, imag);
    std::cout << "time: " << t.elapsed() << "sec" << std::endl;
    for (int i = 0; i< real.size(); i++){
        std::cout << real[i]  << " + i*" << imag[i] << std::endl;
    }
    return 0;
}