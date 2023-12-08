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
    int n=62;

    Mat x = Mat(1, 100);
    Mat y = Mat(1, 100);
    x.scientific(false);
    y.scientific(false);
    //std::vector<double> x;
    //std::vector<double> y;
    Timer t;
    for( int i=1; i<n; i++){
        std::cout << i << std::endl;
        Mat A= makeBigMatrix(i);
        t.reset();
        std::vector<Number> real;
        std::vector<Number> imag;
        eigenvalues_qr_algorithm_givens(A, real, imag);
        //std::cout << i << ": " << t.elapsed() << "sec" << std::endl;
        //y.push_back(t.elapsed());
        //x.push_back(i);
        y[0][i]=t.elapsed();
        x[0][i]=i;
    }

    for (int i=0; i<n; i++){
        std::cout << x[0][i] << ", ";
    }
     std::cout << std::endl;
    for (int i=0; i<n; i++){
        std::cout << y[0][i] << ", " ;
    }
    
    return 0;
}