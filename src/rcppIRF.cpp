// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::field<arma::cube> rcppIRF(arma::cube B, arma::cube H, int R, bool orthogonal = false) {
  int n = B.n_slices;
  int m = B.n_rows;
  int p = B.n_cols/B.n_rows;
  arma::field<arma::cube> IF(n);
  for(int i = 0; i < n; i ++) {
    arma::cube IFi(m, m, R+p, arma::fill::zeros);

    arma::mat I(m, m);
    if(orthogonal) I = chol(H.slice(i));
    else I.eye();

    for(int j = 0; j < p; j++) IFi.slice(j) = I;
    for(int j = p; j < R+p; j++) {
      for(int k = 1; k <= p; k++)
        IFi.slice(j) += B.slice(i).cols((k-1)*m, k*m-1) * IFi.slice(j-k);
    }
    IF(i) = IFi;
  }
  return IF;
}
// arma::field<arma::cube> rcppIRF(arma::cube B, arma::cube H, int R, bool orthogonal = false) {
//   int n = B.n_slices;
//   int m = B.n_cols;
//   arma::field<arma::cube> IF(n);
//   for(int i = 0; i < n; i ++) {
//     arma::cube IFi(m, m, R+1);
//     if(orthogonal) IFi.slice(0) = chol(H.slice(i));
//     else IFi.slice(0).eye();
//     for(int j = 1; j <= R; j++) {
//       IFi.slice(j) = B.slice(i) * IFi.slice(j-1);
//     }
//     IF(i) = IFi;
//   }
//   return IF;
// }
