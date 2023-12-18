#ifndef HDNUM_QR_EIGENVALUES_HH
#define HDNUM_QR_EIGENVALUES_HH

#include "densematrix.hh"
//TODO Eigenvalues < tol set to zero
//TODO value for tol
//TODO description of functions
//TODO tridiagonal matrices (symm matrices in upper hessenberg)

namespace hdnum{

    template <typename T>
    class QR_Info {
        private:
            T tol;
            std::vector<T> real;
            std::vector<T> imag;
        public:
            QR_Info(T tol_): tol(tol_){}

            T get_tol(){
                return tol;
            }
            void add_real(T elem){
                real.push_back(elem);
            }
            void add_imag(T elem){
                imag.push_back(elem);
            }
            std::vector<T> get_real(){
                return real;
            }
            std::vector<T> get_imag(){
                return imag;
            }
    };

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
    void hessenberg_qr(hdnum::DenseMatrix<T>& A, QR_Info<T>& qr_info){
        ///overwrites A with R, old_A=QR; A'=RQ ="AQ"
        int n=A.rowsize();
        T tol=qr_info.get_tol();
        std::vector<std::pair<T, T>> Q;

        for (int j=0; j<n-1; j++){

            if (fabs(A[j+1][j])<=tol){
                A[j+1][j]=0;
                Q.push_back({1, 0});
                continue;
            }
            std::pair<T, T> cs=givens(A[j][j], A[j+1][j]);
            T c=cs.first;
            T s=cs.second;
            Q.push_back(cs);

            for(int k=j; k<n; k++){ //givens^t*A
                T tempAjk=A[j][k];
                A[j][k]=A[j][k]*c-A[j+1][k]*s;
                A[j+1][k]=tempAjk*s+A[j+1][k]*c;
            }
            A[j+1][j]=0; //TODO check this (roundoff correction)
        }

        //multiply with givensroations from the left
        for (int i = 0; i<n-1; i++){  
            for (int k=0; k<n; k++){ //A*givens
                T c=Q[i].first;
                T s=Q[i].second;

                T tempAki=A[k][i];
                A[k][i]=A[k][i]*c-A[k][i+1]*s;
                A[k][i+1]=tempAki*s+A[k][i+1]*c;
            }
        }
    }


    template <class T>
    hdnum::DenseMatrix<T> makeHessenberg(hdnum::DenseMatrix<T>& A, QR_Info<T>& qr_info){ 
        int n=A.rowsize();
        T tol = qr_info.get_tol();
        for(int j=0; j<n-2; j++){
            for(int i=n-2; i>j; i--){
                if (fabs(A[i+1][j])<=tol){
                    A[i+1][j]=0;
                    continue;
                }
                std::pair<T, T> cs=givens(A[i][j], A[i+1][j]);
                T c=cs.first;
                T s=cs.second;

                //A=Qt*A
                for(int k=j; k<n; k++){ //givens^t*A
                    T tempAik=A[i][k];
                    A[i][k]=A[i][k]*c-A[i+1][k]*s;
                    A[i+1][k]=tempAik*s+A[i+1][k]*c;
                }
                //A=A*QT ( --> A=QT*A*Q)
                for (int k=0; k<n; k++){
                    T tempAki=A[k][i];
                    A[k][i]=A[k][i]*c-A[k][i+1]*s;
                    A[k][i+1]=tempAki*s+A[k][i+1]*c;
                }

                A[i+1][j]=0;//roundoff error //TODO check this
            }
        }
        return A;
    }

    template <class T>
    std::vector<int> decouple_check(hdnum::DenseMatrix<T>& A,  QR_Info<T>& qr_info){  //TODO T tol??
        T tol = qr_info.get_tol();
        std::vector<int> zeros;
        for (int i=0; i<A.rowsize()-1; i++){
            if (fabs(A[i+1][i])<tol) {
                A[i+1][i]=0;
                zeros.push_back(i);
            }
        }
        return zeros;
    }

    template <class T>
    void eigenvalues_2x2(const hdnum::DenseMatrix<T>& A,  QR_Info<T>& qr_info){
        assert(A.rowsize()==2 && A.colsize()==2);
        T det=A[0][0]*A[1][1]-A[1][0]*A[0][1];
        T mean=(A[0][0]+A[1][1])/2;
        T root= mean*mean-det;
        if (root>0){
            qr_info.add_real(mean+std::sqrt(root));
            qr_info.add_real(mean-std::sqrt(root));
            qr_info.add_imag(0);
            qr_info.add_imag(0);
            return;
        }
        qr_info.add_real(mean);
        qr_info.add_real(mean);
        qr_info.add_imag(std::sqrt(-1*root));
        qr_info.add_imag(-std::sqrt(-1*root));
        return;
    }

    template <class T>
    void decouple(hdnum::DenseMatrix<T>& A,  QR_Info<T>& qr_info){
        if (A.colsize() < 1) return qr_step(A, qr_info);
        std::vector<int> zeros = decouple_check(A, qr_info);
        if (!zeros.empty()){ //decoupeling
            int n= zeros.size();
            qr_step(A.sub(0, 0, zeros[0]+1, zeros[0]+1), qr_info);
            for (int i=1; i<n; i++){ 
                assert(zeros[i]<A.colsize()-1);
                qr_step(A.sub(zeros[i-1]+1, zeros[i-1]+1, zeros[i]-zeros[i-1], zeros[i]-zeros[i-1]),  qr_info);
            }
            qr_step(A.sub(zeros[n-1]+1, zeros[n-1]+1, A.rowsize()-1-zeros[n-1], A.rowsize()-1-zeros[n-1]), qr_info);
        }
        else qr_step(A, qr_info);
    }

    template <class T>
    void qr_step(hdnum::DenseMatrix<T> A, QR_Info<T>& qr_info){
        if (A.rowsize()==1){
            qr_info.add_real(A[0][0]);
            qr_info.add_imag(0);
        }
        else if (A.rowsize()==2){
            eigenvalues_2x2(A, qr_info);
        }
        else{
            for( int i=0; i<20; i++){ //TODO decouple after every iteration?
                hessenberg_qr(A, qr_info); //= R*Q=A'
            }
            decouple(A, qr_info);
        }
    }


    template <class T> //TODO create real and imag as arrays with size n
    QR_Info<T> eigenvalues_qr_algorithm_givens(hdnum::DenseMatrix<T> A){
        T tol = 1e-16;
        QR_Info<T> qr_info(tol);
        makeHessenberg(A,  qr_info);
        decouple(A, qr_info);
        return qr_info;
    }


} //end namespace hdnum

#endif