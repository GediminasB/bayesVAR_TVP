// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "rcppRmvnorm.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rcppDrawBeta_CC(arma::mat y, arma::cube Z, arma::cube H, arma::cube Q, arma::mat Pi, arma::rowvec beta1, arma::mat P1) {
  // Dimensions
  int n = H.n_cols;
  int m = Q.n_cols;
  int t = H.n_slices;
  // Init matrices
  arma::mat beta_tt(t+1, m);
  arma::cube P_tt(m, m, t+1);
  arma::cube P_t(m, m, t+1);
  arma::rowvec nu(n);
  arma::mat F(n, n); 
  arma::mat Fi(n, n); 
  arma::mat K(m, n); 
  arma::mat L(m, m); 
  // Prior for a1
  beta_tt.row(0) = beta1;
  P_t.slice(0) = P1;
  // Kalman filter
  for(int i = 0; i < t; i++) {
    nu = y.row(i) - beta_tt.row(i) * Z.slice(i).t();
    F = Z.slice(i) * P_t.slice(i) * Z.slice(i).t() + H.slice(i);
    Fi = F.i();
    K = Pi * P_t.slice(i) * Z.slice(i).t() * Fi;
    L = Pi - K * Z.slice(i);
    beta_tt.row(i+1) = beta_tt.row(i) * Pi.t() + nu * K.t();
    P_tt.slice(i+1) = Pi * P_t.slice(i) * L.t();
    P_t.slice(i+1) = P_tt.slice(i+1) + Q.slice(i);
  }
  
  arma::mat beta_draw(t, m, arma::fill::zeros);
  beta_draw.row(t-1) = rcppRmvnorm(1, beta_tt.row(t).t(), P_tt.slice(t));   
  // Backward recursions
  for (int i = t-1; i > 0; i--) {
    arma::mat P_ti = P_t.slice(i).i();
    arma::rowvec beta_mean = beta_tt.row(i) + (beta_draw.row(i) - beta_tt.row(i)) * (P_tt.slice(i)*P_ti).t();
    arma::mat beta_var = P_tt.slice(i) - P_tt.slice(i) * P_ti * P_tt.slice(i);
    beta_draw.row(i-1) = rcppRmvnorm(1, beta_mean.t(), beta_var);
  }
  return beta_draw;
}