#ifndef HDNUM_QR_EIGENVALUES_HH
#define HDNUM_QR_EIGENVALUES_HH

#include "densematrix.hh"
//TODO Eigenvalues < tol set to zero
//TODO value for tol
//TODO description of functions
//TODO use efficent algorithm for R*Q (R upper triangle, Q upper hessenberg) --> + Q=G1*G2*... inefficient ?
//TODO tridiagonal matrices (symm matrices in upper hessenberg)

namespace hdnum{
    template <class T> //TODO change insert_partial_matrix to fit parameters of sub
    void insert_partial_matrix(hdnum::DenseMatrix<T>& A,const hdnum::DenseMatrix<T>& Partial, int row_start, int row_end, int col_start, int col_end){
        assert(row_end<A.rowsize() && col_end<A.colsize());
        assert(row_end-row_start == Partial.rowsize()-1 && col_end-col_start == Partial.colsize()-1);
        for (int i=col_start; i<=col_end; i++){
            for (int j=row_start; j<=row_end; j++){
                A[j][i] = Partial[j-row_start][i-col_start];
            }
        }
    }

    template <class T>
    std::pair<T, T> givens(T a, T b){
        if(b==0){
            return std::pair<T, T> (1, 0);
        }
        T t;
        T c;
        T s;
        if(abs(b)>abs(a)){
            t = -a/b;
            s = 1/(std::sqrt(1+t*t));
            c = s*t;
        }
        else{
            t = -b/a;
            c = 1/(std::sqrt(1+t*t));
            s = c*t;
        }
        return std::pair<T, T> (c, s);
    }

    template <class T>
    void hessenberg_qr(hdnum::DenseMatrix<T>& A, T tol){
        ///overwrites A with R, old_A=QR; A'=RQ ="AQ"
        int n=A.rowsize();
        hdnum::DenseMatrix<T> Q(n, n);
        for(int i =0; i<n; i++) Q[i][i]=1;

        for (int j=0; j<n-1; j++){
            if (fabs(A[j+1][j])<=tol){
                A[j+1][j]=0;
                continue;
            }
            std::pair<T, T> cs=givens(A[j][j], A[j+1][j]);
            T c=cs.first;
            T s=cs.second;
            hdnum::DenseMatrix<T> givens_part={{c, s}, {-s, c}}; //TODO Segmentation faul for some values (e.g. 92)
            //use multiplication facor instead of matrix
            hdnum::DenseMatrix<T> a_part=A.sub(j, j, 2, n-j);
            a_part=givens_part.transpose()*a_part;
            insert_partial_matrix(A, a_part, j, j+1, j, n-1);
            A[j+1][j]=0; //TODO check this (roundoff correction)

            hdnum::DenseMatrix<T> q_part=Q.sub(0, j, j+2, 2);
            q_part=q_part*givens_part;
            insert_partial_matrix(Q, q_part, 0, j+1, j, j+1);
        }
        A=A*Q;
    }


    template <class T>
    hdnum::DenseMatrix<T> makeHessenberg(hdnum::DenseMatrix<T>& A, T tol){ 
        int n=A.rowsize();
        for(int j=0; j<n-2; j++){
            for(int i=n-2; i>j; i--){
                if (fabs(A[i+1][j])<=tol){
                    A[j+1][j]=0;
                    continue;
                }
                std::pair<T, T> cs=givens(A[i][j], A[i+1][j]);
                T c=cs.first;
                T s=cs.second;
                hdnum::DenseMatrix<T> givens_part={{c, s}, {-s, c}};

                //A=Qt*A
                hdnum::DenseMatrix<T> a_part=A.sub(i, j, 2, n-j); //rows
                a_part=givens_part.transpose()*a_part;
                insert_partial_matrix(A, a_part, i, i+1, j, n-1);

                //A=A*QT ( --> A=QT*A*Q)
                a_part=A.sub(0, i, n, 2); //columns
                a_part=a_part*givens_part;
                insert_partial_matrix(A, a_part, 0, n-1, i, i+1);
                //std::cout << "A" << i+1 << ", " << i << " = " << A[i+1][j] << std::endl; 
                A[i+1][j]=0;//roundoff error //TODO check this
            }
        }
        return A;
    }

    template <class T>
    std::vector<int> decouple_check(hdnum::DenseMatrix<T> A, T tol){  //TODO T tol??
        std::vector<int> zeros;
        for (int i=0; i<A.rowsize()-1; i++){
            if (fabs(A[i+1][i])<tol) {zeros.push_back(i);};
        }
        return zeros;
    }

    template <class T>
    void eigenvalues_2x2(const hdnum::DenseMatrix<T>& A, std::vector<T>& real, std::vector<T>& imag){
        T det=A[0][0]*A[1][1]-A[1][0]*A[0][1];
        T mean=(A[0][0]+A[1][1])/2;
        T root= mean*mean-det;
        if (root>0){
            real.push_back(mean+std::sqrt(root));
            real.push_back(mean-std::sqrt(root));
            imag.push_back(0);
            imag.push_back(0);
            return;
        }
        real.push_back(mean);
        real.push_back(mean);
        imag.push_back(std::sqrt(-1*root));
        imag.push_back(-std::sqrt(-1*root));
        return;
    }

    template <class T>
    void decouble(hdnum::DenseMatrix<T> A, T tol,  std::vector<T>& real, std::vector<T>& imag){
        std::vector<int> zeros = decouple_check(A, tol);
        if (!zeros.empty()){ //decoubeling
            int n= zeros.size();
            if(n>0){
                qr_iteration(A.sub(0, 0, zeros[0]+1, zeros[0]+1), tol, false, real, imag);
            }
            for (int i=1; i<n; i++){ 
                qr_iteration(A.sub(zeros[i-1]+1, zeros[i-1]+1, zeros[i]-zeros[i-1], zeros[i]-zeros[i-1]), tol, false, real, imag);
            }
            qr_iteration(A.sub(zeros[n-1]+1, zeros[n-1]+1, A.rowsize()-1-zeros[n-1], A.rowsize()-1-zeros[n-1]), tol, false, real, imag);
        }
        else qr_step(A, tol, real, imag);
    }

    template <class T>
    void qr_step(hdnum::DenseMatrix<T> A, T tol, std::vector<T>& real, std::vector<T>& imag){
        if (A.rowsize()==1){
            real.push_back(A[0][0]);
            imag.push_back(0);
        }
        else if (A.rowsize()==2){
            eigenvalues_2x2(A, real, imag);
        }
        else{
            for( int i=0; i<20; i++){ //TODO decouple after every iteration?
                hessenberg_qr(A, tol); //= R*Q=A'
            }
            qr_iteration(A, tol, true, real, imag);
        }
    }

    template <class T>
    void qr_iteration(hdnum::DenseMatrix<T> A, T tol, bool decouple_bool, std::vector<T>& real, std::vector<T>& imag){
        if (decouple_bool){    
            decouble(A, tol, real, imag);
            }
        else qr_step(A, tol, real, imag);
    }

    template <class T> //TODO create real and imag as arrays with size n
    void eigenvalues_qr_algorithm_givens(hdnum::DenseMatrix<T> A, std::vector<T>& real, std::vector<T>& imag ){
        T tol = 1e-16;
        makeHessenberg(A, tol);
        qr_iteration(A, tol, true, real, imag);
    }



} //end namespace hdnum

#endif