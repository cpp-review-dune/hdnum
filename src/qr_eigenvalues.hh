#ifndef HDNUM_QR_EIGENVALUES_HH
#define HDNUM_QR_EIGENVALUES_HH

#include "densematrix.hh"


namespace hdnum{
    template <class T>
    void insert_partial_matrix(hdnum::DenseMatrix<T>& A,hdnum::DenseMatrix<T> Partial, int row_start, int row_end, int col_start, int col_end){
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
    hdnum::DenseMatrix<T> hessenberg_qr(hdnum::DenseMatrix<T>& A){
        ///overwrites A with R, old_A=QR; A'=RQ =AQ; returns Q
        int n=A.rowsize();
        hdnum::DenseMatrix<T> Q(n, n);
        for(int i =0; i<n; i++) Q[i][i]=1;

        for (int j=0; j<n-1; j++){
            std::pair<T, T> cs=givens(A[j][j], A[j+1][j]);
            T c=cs.first;
            T s=cs.second;
            hdnum::DenseMatrix<T> givens_part={{c, s}, {-s, c}};

            hdnum::DenseMatrix<T> a_part=A.sub(j, j, 2, n-j);
            a_part=givens_part.transpose()*a_part;
            insert_partial_matrix(A, a_part, j, j+1, j, n-1);

            hdnum::DenseMatrix<T> q_part=Q.sub(0, j, n, 2);
            q_part=q_part*givens_part;
            insert_partial_matrix(Q, q_part, 0, n-1, j, j+1);
        }
        return Q;
    }


    template <class T>
    hdnum::DenseMatrix<T> makeHessenberg(hdnum::DenseMatrix<T>& A){ 
        int n=A.rowsize();
        for(int j=0; j<n-2; j++){
            for(int i=n-2; i>j; i--){
                if (abs(A[i+1][j]) <= 0.00000000000001) continue; //TODO add tol
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
            }
        }
        return A;
    }

    template <class T>
    std::vector<int> decouple_check(hdnum::DenseMatrix<T> A, T tol){  //TODO T tol??
        std::vector<int> zeros({});
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
    void qr_iteration(hdnum::DenseMatrix<T> A, T tol, bool decouple, std::vector<T>& real, std::vector<T>& imag){
        bool qr_step=false;
        if (decouple){    
            std::vector<int> zeros = decouple_check(A, tol);
            if (!zeros.empty()){ //decoubeling
                int n= zeros.size();
                for (int i=0; i<n; i++){
                    if(i==0) qr_iteration(A.sub(0, 0, zeros[i]+1, zeros[i]+1), tol, false, real, imag);
                    else{
                        qr_iteration(A.sub(zeros[i-1]+1, zeros[i-1]+1, zeros[i]-zeros[i-1], zeros[i]-zeros[i-1]), tol, false, real, imag);
                    }
                }
                qr_iteration(A.sub(zeros[n-1]+1, zeros[n-1]+1, A.rowsize()-1-zeros[n-1], A.rowsize()-1-zeros[n-1]), tol, false, real, imag);

            }
            else qr_step = true;
        }
        else qr_step = true;
        if (qr_step){
            if (A.rowsize()==1){
                real.push_back(A[0][0]);
                imag.push_back(0);
            }
            else if (A.rowsize()==2){
                eigenvalues_2x2(A, real, imag);
            }
            else{
                for( int i=0; i<20; i++){ //TODO decouple after every iteration?
                    A=A* hessenberg_qr(A); //= R*Q=A'
                }
                qr_iteration(A, tol, true, real, imag);
            }
        }
    }

    template <class T>
    void eigenvalues_qr_algorithm_givens(hdnum::DenseMatrix<T> A, std::vector<T>& real, std::vector<T>& imag ){
        T tol = 0.0000000000000001;
        makeHessenberg(A);
        qr_iteration(A, tol, true, real, imag);
    }



} //end namespace hdnum

#endif