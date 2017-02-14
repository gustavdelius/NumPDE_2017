#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector doublesweep(NumericVector A, NumericVector B, NumericVector C, 
                          NumericVector F, double a, double b) {
    int N = C.size() + 1;
    if ((B.size() != N-1) || (A.size() != N-1) || (F.size() != N-1)) {
        stop("The lengths of the vector arguments need to be equal");
    }
    
    NumericVector alpha(N);
    NumericVector beta(N);
    alpha[0] = 0;
    beta[0] = a;
        
    for (int i=0; i < (N-1); ++i) {
        alpha[i+1] = B[i] / (C[i]-alpha[i]*A[i]);
        beta[i+1] = (beta[i]*A[i] - F[i]) / (C[i] - alpha[i]*A[i]);
    }
    
    NumericVector v(N-1);
    v[N-2] = alpha[N-1]*b + beta[N-1];
    
    for (int i=N-2; i > 0; --i) {
        v[i-1] = alpha[i]*v[i] + beta[i];
    }

    return v;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
A <- c(0,2,3,1)
B <- c(2,1,4,0)
C <- c(1,1,1,7)
F <- c(1,2,6,-6)
doublesweep(A,B,C,F,0,0)
*/
