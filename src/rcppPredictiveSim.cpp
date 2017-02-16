// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "rcppRmvnorm.h"

// [[Rcpp::export]]
arma::cube rcppMakeRegs(arma::mat y, int p = 1) {
  int T = y.n_rows;
  int n = y.n_cols;
  arma::cube Z(n, n*(n*p+1), T-p);
  arma::mat I(n, n, arma::fill::eye);
  arma::vec one(1, arma::fill::ones);
  arma::mat yt = fliplr(y.t());
  for(int i = 0; i < T-p; i++) {
    Z.slice(i) = arma::kron(arma::join_cols(one, arma::vectorise(yt.cols(T-i-p, T-i-1))).t(), I);
  }
  return(Z);
}
arma::mat rcppMakeRegs_one(arma::mat y) {
  int n = y.n_cols;
  arma::mat I(n, n, arma::fill::eye);
  arma::vec one(1, arma::fill::ones);
  arma::mat yt = fliplr(y.t());
  return(arma::kron(arma::join_cols(one, arma::vectorise(yt)).t(), I));
}
// [[Rcpp::export]]
arma::cube rcppPredictiveSim_TVP(arma::mat starting_Y, arma::mat B, arma::cube H, arma::cube Q, int T) {
  int N = B.n_rows;
  int n_vars = B.n_cols;
  int n = H.n_cols;
  int p = (n_vars - n)/n/n;
  arma::cube F(T+p, n, N);

  for(int i = 0; i < N; i++) {
    arma::rowvec beta = B.row(i);
    arma::mat yF(T, n);
    yF.insert_rows(0, starting_Y);
    //yF.rows(0, p-1) = starting_Y;
    for(int j = p; j < T+p; j++) {
      beta = beta + rcppRmvnorm(1, arma::zeros(n_vars), Q.slice(i));
      yF.row(j) = beta * rcppMakeRegs_one(yF.rows(j-p, j-1)).t() + rcppRmvnorm(1, arma::zeros(n), H.slice(i));
    }
    F.slice(i) = yF;
  }
  return F;
}
// [[Rcpp::export]]
arma::cube rcppPredictiveSim_FTVP(arma::mat starting_Y, arma::mat B, arma::cube H, arma::cube Q, arma::cube L, arma::cube R, int T, int k) {
  int N = B.n_rows;
  int n_vars = B.n_cols;
  int n = H.n_cols;
  int p = (n_vars - n)/n/n;
  int l = L.n_cols;
  int n_f = L.n_rows;
  int n_y = starting_Y.n_cols;
  arma::cube F(T+p, n_y, N*k);

  for(int i = 0; i < N; i++) {
    arma::mat S = L.slice(i).t() * (L.slice(i) * L.slice(i).t()).i();
    arma::mat starting_F = join_rows(starting_Y.cols(0, l-1) * S, starting_Y.cols(l, n_y-1));
    for(int ii = 0; ii < k; ii++) {
      arma::rowvec beta = B.row(i);
      arma::mat FF(T, n);
      arma::mat yF(T, n_y);
      FF.insert_rows(0, starting_F);
      yF.insert_rows(0, starting_Y);
      //yF.rows(0, p-1) = starting_Y;
      for(int j = p; j < T+p; j++) {
        beta = beta + rcppRmvnorm(1, arma::zeros(n_vars), Q.slice(i));
        FF.row(j) = beta * rcppMakeRegs_one(FF.rows(j-p, j-1)).t() + rcppRmvnorm(1, arma::zeros(n), H.slice(i));
        yF.row(j) = join_rows(FF.row(j).cols(0,n_f-1) * L.slice(i) + rcppRmvnorm(1, arma::zeros(R.n_cols), R.slice(i)), FF.row(j).cols(n_f,n-1));
      }
      F.slice(k*i + ii) = yF;
    }
  }
  return F;
}
