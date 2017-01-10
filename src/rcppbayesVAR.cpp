// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::cube rcppExpandKronecker(arma::mat x, int n) {
  int t = x.n_rows;
  int m = x.n_cols * n;
  arma::cube Z(n, m, t);
  arma::mat I(n, n, arma::fill::eye);
  for(int i = 0; i < t; i++) {
    Z.slice(i) = arma::kron(x.row(i), I);
  }
  return Z;
}

// [[Rcpp::export]]
arma::mat rcppZHZ(arma::cube Z, arma::mat H) {
  int n = Z.n_slices;
  int m = Z.n_cols;
  arma::mat ZHZ(m, m, arma::fill::zeros);
  for(int i = 0; i < n; ++i) {
    ZHZ = ZHZ + Z.slice(i).t() * H * Z.slice(i);
  }
  return ZHZ;
}

// [[Rcpp::export]]
arma::mat rcppZHy(arma::cube Z, arma::mat H, arma::mat y) {
  int n = Z.n_slices;
  int m = Z.n_cols;
  arma::mat ZHy(m, 1, arma::fill::zeros);
  for(int i = 0; i < n; ++i) {
    ZHy = ZHy + Z.slice(i).t() * H * y.row(i).t();
  }
  return ZHy;
}

// [[Rcpp::export]]
arma::mat rcppSSE(arma::mat y, arma::cube Z, arma::colvec beta) {
  int n = Z.n_slices;
  int m = y.n_cols;
  arma::mat SSE(m, m, arma::fill::zeros);
  for(int i = 0; i < n; ++i) {
    SSE = SSE + (y.row(i).t() - Z.slice(i) * beta) * (y.row(i).t() - Z.slice(i) * beta).t();
  }
  return SSE;
}

// [[Rcpp::export]]
arma::mat rcppSSEmat(arma::mat y, arma::cube Z, arma::mat beta) {
  int n = Z.n_slices;
  int m = y.n_cols;
  arma::mat SSE(m, m, arma::fill::zeros);
  for(int i = 0; i < n; ++i) {
    SSE = SSE + (y.row(i).t() - Z.slice(i) * beta.row(i).t()) * (y.row(i).t() - Z.slice(i) * beta.row(i).t()).t();
  }
  return SSE;
}

// [[Rcpp::export]]
arma::mat rcppResidConstB(arma::mat y, arma::cube Z, arma::rowvec beta) {
  int n = Z.n_slices;
  arma::mat Resid(y.n_rows, y.n_cols);
  for(int i = 0; i < n; ++i) {
    Resid.row(i) = y.row(i) - beta * Z.slice(i).t();
  }
  
  return Resid;
}
// [[Rcpp::export]]
arma::field<arma::cube> rcppIRF(arma::cube B, arma::cube H, int R, bool orthogonal = false) {
  int n = B.n_slices;
  int m = B.n_cols;
  arma::field<arma::cube> IF(n);
  for(int i = 0; i < n; i ++) {
    arma::cube IFi(m, m, R+1);
    if(orthogonal) IFi.slice(0) = chol(H.slice(i));
    else IFi.slice(0).eye();
    for(int j = 1; j <= R; j++) {
      IFi.slice(j) = B.slice(i) * IFi.slice(j-1); 
    }
    IF(i) = IFi;
  }
  return IF;
}