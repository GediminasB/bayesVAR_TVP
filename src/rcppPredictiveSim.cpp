// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "rcppRmvnorm.h"

// [[Rcpp::export]]
arma::cube rcppPredictiveSim_TVP(arma::rowvec starting_Y, arma::mat B, arma::cube H, arma::cube Q, int T) {
  int N = B.n_rows;
  int n_vars = B.n_cols;
  int n = H.n_cols;
  arma::cube F(T+1, n, N);
  arma::mat I(n, n, arma::fill::eye);

  for(int i = 0; i < N; i++) {
    arma::rowvec beta = B.row(i);
    arma::mat yF(T+1, n);
    yF.row(0) = starting_Y;
    for(int j = 1; j <= T; j++) {
      beta = beta + rcppRmvnorm(1, arma::zeros(n_vars), Q.slice(i));
      yF.row(j) = beta.subvec(0, n-1) + beta.subvec(n, n_vars-1) * arma::kron(yF.row(j-1), I).t() + rcppRmvnorm(1, arma::zeros(n), H.slice(i));
    }
    F.slice(i) = yF;
  }
  return F;
}
