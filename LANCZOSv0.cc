#include <iostream>
#include <vector>
#include <valarray>
//#include "hdnum.hh"

//TODO: Einbindung hdnum
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

void lanzcos(std::vector<std::vector<double>> A, std::vector<double> q1){
    /**
     * initialisation:
     * k = 0; \beta_0 = 1, q_0 = 0, r_0 = q_1
    */

    int k = 0;
    std::vector<double> beta = {1}; //beta0 = 1
    std::vector<double> q0;
    for (int i = 0; i < A.size(); i++) {
        q0.push_back(0);
    }
    std::vector<std::vector<double>> Q = {q0};
    std::vector<double> r = q1; 
    std::vector<double> alpha;
    
    while (k == 0 || beta[k] != 0) {
        // q_(k+1) = r_k / \beta_k
        std::vector<double> q_kp1 = r;
        for (auto i : q_kp1) {
            i = i/beta[k];
        }
        Q.push_back(q_kp1);

        k++;
        // alpha_k = q_k^T * A * q_k
        double alpha_k; //TODO: alpha_k = q_k^T * A * q_k (matrix multiplication)

        // r = (A - alpha_k * I) * q[k] - beta[k-1] * q[k-1]
        r; //TODO: r = (A - alpha_k * I) * q[k] - beta[k-1] * q[k-1] (matrix multiplication)

        // beta[k] = 2-norm of r
        // var1
        double underroot = 0;
        for (auto i : r) {
            underroot = underroot + std::pow(std::abs(i), 2);
        }
        //var2: transpose(r) * r
        underroot = 0;
        for (auto i : r) {
            underroot = underroot + i * i;
        }
        beta[k] = std::sqrt(underroot);
    }

    //TODO: return matrix T_k
}

int main(){
    //lanzcos();

    return 0;
}