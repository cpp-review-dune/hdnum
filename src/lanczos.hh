// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef LANCZOS_HH
#define LANCZOS_HH

#include "vector.hh"
#include "densematrix.hh"
#include "sparsematrix.hh"
#include <random>
#include <cfloat>

//TODO: termination condition: 
    //TODO: set precision to more than 3
    //USE: std::abs(beta[k]) > 1e-9
//TODO: eigenvalues von Tk -> ausnutzen Bandmatrix
//TODO: test application without randomized q_1

///UPDATE: storage optimization
///UPDATE: sparseMatrix included
///UPDATE: first CRO method implemented

/** @file
 *  @brief This file implements the Lanczos method
 */

/**
 * input: sparse symmetric matrix A \in \R^(n \times n), unit 2-norm vector q_1 \in \R^n 
 * output: matrix T_k \in R^(k \times k)
 * 
 * Data structures:
 * - matrix A \in \R^(n \times n)
 * - matrix Q \in \R^(n \times k) consists of vectors q_i \in \R^n, 1 \le i \le k 
 * - scalars \alpha_i, 1 \le i \le k
 * - scalars \beta_i, 1 \le i \le k
 * - integer k, 1 \le k \le n
*/
namespace hdnum {

/************************************************************************************************/
/* UTILITIES */
/************************************************************************************************/

    //Number comparison
    //code from : https://stackoverflow.com/questions/19837576/comparing-floating-point-number-to-zero
    template<typename T>
    bool approx_equal(const T &a, const T &b){
        if (a == b) return true; 
        return false;
    }

    template<>
    bool approx_equal<double>(const double &a, const double &b){
        if (std::fabs(a - b) <= DBL_EPSILON * std::fmax(std::fabs(a), std::fabs(b))) return true; 
        return false;
    }

    template<>
    bool approx_equal<long double>(const long double &a, const long double &b){
        if (std::fabs(a - b) <= LDBL_EPSILON * std::fmax(std::fabs(a), std::fabs(b))) return true; 
        return false;
    }
    template<>
    bool approx_equal<float>(const float &a, const float &b){
        if (std::fabs(a - b) <= FLT_EPSILON * std::fmax(std::fabs(a), std::fabs(b))) return true; 
        return false;
    }

    template<class T>
    bool approx_a_less_b(const T &a, const T &b){
        T amb = a - b;
        if (amb < 0 && !approx_equal(amb, (T) 0)) return true;
        return false;
    }

    template<class T>
    bool approx_a_lesseq_b(const T &a, const T &b){
        T amb = a - b;
        if (amb < 0 || approx_equal(amb, (T) 0)) return true;
        return false;
    }

    
    //TODO: make v a matrix instead of a Vector 
    template<class T>
    std::pair<Vector<T>, T> householder(Vector<T> x) {
        T beta;
        //m = length(x); sigma = x(2:m)^T * x(2:m); v = [1, x(2:m)]
        size_t m = x.size();
        Vector<T> x_sub(m-1); x_sub = x.sub(1, m-1);
        std::cout << "x_sub: " << x_sub << std::endl;
        T sigma = x_sub * x_sub;
        Vector<T> v(m); v[0] = 1;
        for (size_t i = 1; i < m; i++) {
            v[i] = x[i];
        }
        std::cout << sigma << std::endl;
        //if (sigma = 0 && x(1) >= 0) beta = 0
        //else if (sigma = 0 && x(1) < 0) beta = -2
        //TODO: test for zero
        if (approx_equal(sigma, 0.0)){
            //x[0] < 0 
            if (approx_a_less_b(x[0], 0.0)) {
                beta = -2;
            }
            else {
                beta = 0;
            }
        }
        else {
            //mu = sqrt(x(1)**2 + sigma)
            T mu = std::sqrt(std::pow(x[0], 2) + sigma);
            //if (x(1) <= 0) v(1) = x(1) - mu
            //else v(1) = - sigma/(x(1) + mu)
            v[0] = sigma/(x[0] + mu);
            std::cout << v << std::endl;
            if (x[0] > 0) {
                v[0] *= -1;
            }
            //beta = 2 * v(1)**2 / (sigma + v(1) ** 2)
            beta = 2 * std::pow(v[0], 2) / (sigma + std::pow(v[0], 2));
            //v = v/v(1)
            v /= v[0];
        }
        std::cout << beta << std::endl;
        return std::make_pair(v, beta);
    }

    /**
     * @brief random q_1 vector generator
     * 
     * @param n  size of matrix \f$ A \f$
     * @return Vector<T> 2-norm unit vector \f$ \in \R^n \f$
     */
    //TODO: make template -> pt everything into class 
    Vector<double> generate_q_1(const size_t &n) {
        //random double generated according to https://www.geeksforgeeks.org/generate-random-T-numbers-in-cpp/
        double lower_bound = -1;
        double upper_bound = 1;
        std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
        std::default_random_engine re;
        
        //init random vector
        Vector<double> q_1 (n);
        for (size_t i = 0; i < n; i++) {
            q_1[i] = unif(re);
        }

        //normalize to 2-norm
        q_1 /= q_1.two_norm();

        return q_1;
    }

/************************************************************************************************/
/* LANCZOS BASIC */
/************************************************************************************************/

    /**
     * @brief basic algorithm according to book with no storage or numerical optimization
     * 
     * @tparam T 
     * @param A 
     * @return SparseMatrix<T> 
     */
    template<class T>
    SparseMatrix<T> lanczos_basic(const SparseMatrix<T> &A, const Vector<T> &q_1){
        //init:
        int n = A.rowsize();
        int k = 0;
        std::vector<T> beta {1}; //beta_0 = 1
        Vector<T> q_km1(n, 0); //q_0 = 0
        Vector<T> r(n);
        r = q_1; //r_0 = generate_q_1
        std::vector<T> alpha;

        //TODO: additional termination condition
        while (!approx_equal(beta[k], 0.0) && k < n) {
            k++;
            // q_(k) = r_(k-1) / \beta_(k-1)
            auto q_k = r;
            q_k /= beta[k-1];

            // alpha_k = q_k^T * A * q_k
            Vector<T> Aq_k(n);
            A.mv(Aq_k, q_k);
            auto alpha_k = q_k * Aq_k;
            alpha.push_back(alpha_k);

            // r = (A - alpha_k * I) * q[k] - beta[k-1] * q[k-1]
            //(A - alpha_k * I) shortened to A(i,i) - alpha_k
            auto Amalpha = A;
            for (int i = 0; i < Amalpha.rowsize(); i++) {
                Amalpha.get(i,i) = Amalpha(i,i) - alpha_k;
            }
            Amalpha.mv(r, q_k);
            auto betaq = q_km1;
            betaq *= beta[k-1];
            r -= betaq; 

            // beta[k] = 2-norm of r
            beta.push_back(r.two_norm());
            q_km1 = q_k; 
        }

        SparseMatrix<T> Tk (k,k);
        auto builder = typename SparseMatrix<T>::builder(k,k);
        for (size_t i = 0; i < k; i++) {
            builder.addEntry(i, i, alpha[i]);
            if (i != (k-1)) {
                builder.addEntry(i, i+1, beta[i]);
                builder.addEntry(i+1, i, beta[i]);
            }
        }
        Tk = builder.build();
        return Tk;
    }

    template<class T>
    SparseMatrix<T> lanczos_basic(const SparseMatrix<T> &A){
        SparseMatrix<T> Tk;
        Tk = lanczos_basic(A, generate_q_1(A.colsize())); 
        return Tk;
    }

/************************************************************************************************/
/* LANCZOS OPTIMIZED */
/************************************************************************************************/

    /**
     * @brief computes a matrix \f$ T_k \f$ thats extremal eigenvalues are approximations 
     *        for the extreme eigenvalues of a large, sparse, symmetric matrix \f$ A \f$
     * 
     * @tparam T 
     * @param A a (large,) sparse, symmetric matrix 
     * @return SparseMatrix<T> matrix T_k
     */
    template<class T>
    SparseMatrix<T> lanczos_optimized(const SparseMatrix<T> &A, const Vector<T> &q_1) {
        //init
        int n = A.rowsize();
        //k = 1 equals k = 0 because our numeration starts at 0
        int k = 0;

        Vector<T> w(n);
        Vector<T> v(n);
        std::vector<T> alpha;
        std::vector<T> beta;

        //first iteration:
        //w = q_1 
        w = q_1;
        //v = A * w
        A.mv(v, w);
        //\alpha_1 = w^T * v
        alpha.push_back(w * v);
        //v = v - \alpha_1 * w
        Vector<T> alphaw(n);
        alphaw = w;
        alphaw *= alpha[0];
        v -= alphaw;
        //\beta_1 = ||v||_2
        beta.push_back(v.two_norm());

        // k <= 2*n for termination
        //TODO: termination
        while (!approx_equal(beta[k], 0.0) && k < n) {
            for (int i = 0; i < n; i++) {
                //t = w_i; w_i = v_i/\beta_k; v_i = -\beta_k * t
                T t = w[i];
                w[i] = v[i] / beta[k];
                v[i] = -1 * beta[k] * t;
            }

            //v = v + A*w
            Vector<T> Aw(n);
            A.mv(Aw, w);
            v += Aw;
            //
            k++;
            //\alpha_k = w^T * v
            alpha.push_back(w * v);
            //v = v - \alpha_k * w
            alphaw = w;
            alphaw *= alpha[0];
            v -= alphaw;
            //\beta_k = ||v||_2
            beta.push_back(v.two_norm());
        }

        // DenseMatrix<T> Tk(k, k);
        // for(int i = 0; i < k; i++) {
        //     Tk[i][i] = alpha[i];
        //     if (i != (k-1)) {
        //         Tk[i][i+1] = beta[i];
        //         Tk[i+1][i] = beta[i];
        //     }
        // }

        SparseMatrix<T> Tk;
        auto builder = typename SparseMatrix<T>::builder(k,k);
        for (size_t i = 0; i < k; i++) {
            builder.addEntry(i, i, alpha[i]);
            if (i != (k-1)) {
                builder.addEntry(i, i+1, beta[i]);
                builder.addEntry(i+1, i, beta[i]);
            }
        }
        Tk = builder.build();
        return Tk;
    }

    template<class T>
    SparseMatrix<T> lanczos_optimized(const SparseMatrix<T> &A){
        SparseMatrix<T> Tk;
        Tk = lanczos_optimized(A, generate_q_1(A.colsize())); 
        return Tk;
    }

/************************************************************************************************/
/* LANCZOS COMPLETE REORTHOGONALIZATION */
/************************************************************************************************/

    /**
     * @brief OTHER ALGORITHM ATTEMPT
     * 
     * @tparam T 
     * @param A 
     * @return SparseMatrix<T> 
     */
    template<class T>
    SparseMatrix<T> lanczos_cro(const SparseMatrix<T> &A, const Vector<T> &q_1) {
        //init:
        assert(A.rowsize() == A.colsize());

        int n = A.rowsize();
        int k = 0;

        std::vector<T> alpha;
        std::vector<T> beta;
        
        Vector<T> q(n);
        DenseMatrix<T> Q;
        DenseMatrix<T> QT;
        Vector<T> r(n);

        //q = x/||x||; Q = [q]
        q = q_1;
        QT = Q.transpose(); QT.addNewRow(q_1); Q = QT.transpose();

        //first iteration
        //r = A * q
        A.mv(r, q);
        //alpha_1 = q * r
        alpha.push_back(q * r);
        //r = r - alpha_1 * q
        r.update(-1 * alpha.back(), q);
        //beta_1 = ||r||
        beta.push_back(r.two_norm());
        T beta_k = r.two_norm();
        
        while(fabs(beta[k]) > 1e-18 && k < n) {
            k++;
            //v = q_(k-1)
            auto q_km1 = q;

            //q = r/beta_(j-1)
            q = r;
            q /= beta[k-1];

            //Q_j = [Q_(j-1), q]
            QT = Q.transpose(); QT.addNewRow(q); Q = QT.transpose();

            //r = A * q - beta_(j-1) * v
            A.mv(r, q);
            r.update(-1 * beta.back(), q_km1);

            //alpha_j = q * r
            alpha.push_back(q * r);

            //r = r - alpha_j * q
            r.update(-1 * alpha.back(), q);

            //reorthogonalization: r = r - Q * (Q^T * r)
            Vector<T> QTr (Q.colsize());
            QT.mv(QTr, r);
            Vector<T> QQTr (n);
            Q.mv(QQTr, QTr);
            r -= QQTr;

            //beta_j = ||r||
            beta_k = r.two_norm();
            beta.push_back(r.two_norm());
        }

        std::cout << k << std::endl;

        SparseMatrix<T> Tk;
        auto builder = typename SparseMatrix<T>::builder(k,k);
        for (int i = 0; i < k; i++) {
            builder.addEntry(i, i, alpha[i]);
            if (i != (k-1)) {
                builder.addEntry(i, i+1, beta[i]);
                builder.addEntry(i+1, i, beta[i]);
            }
        }
        Tk = builder.build();

        return Tk;
    }

    template<class T>
    SparseMatrix<T> lanczos_cro(const SparseMatrix<T> &A){
        SparseMatrix<T> Tk;
        Tk = lanczos_cro(A, generate_q_1(A.colsize())); 
        return Tk;
    }

/************************************************************************************************/
/* LANCZOS COMPLETE REORTHOGONALIZATION USING HOUSEHOLDER */
/************************************************************************************************/

    /**
     * @brief Lanczos method with complete reorthogonalization using Householder matrices
     * 
     * @tparam T 
     * @param A 
     * @return SparseMatrix<T> 
     */
    // template<class T>
    // SparseMatrix<T> lanczos_householder(const SparseMatrix<T> &A, const Vector<T> &q_1) {
    //     /*
    //      * Algorithm:
    //      * 
    //      * r_0 = q_1 (given unit vector)
    //      * Determine Householder H_0 s.t. H_0 * r_0 = e_1
    //      * for (k = 1 .. n-1) {
    //      *      a_k = q_k^T * A * q_k
    //      *      r_k = (A-alpha_k * I) * q_k - beta_(k-1) * q_(k-1)
    //      *      w = (H_(k-1) * ... * H_0) * r_k
    //      *      Determine Householder H_k s.t. H_k * w = [w_1, ..., w_k, beta_k, 0,...,0]^T
    //      *      q_(k+1) = H_0 * .... * H_k * e_(k+1)
    //      * }
    //      */


    //     SparseMatrix<T> Tk(k,k);
    //     return T_k;
    // }

/************************************************************************************************/
/* LANCZOS SELECTIVE REORTHOGONALIZATION */
/************************************************************************************************/


} //namespace hdnum
#endif //LANCZOS_HH