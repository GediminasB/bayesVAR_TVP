// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "rcppRmvnorm.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rcppKalmanFilterSmooth(arma::mat y, arma::cube Z, arma::cube H, arma::cube Q, arma::mat Pi, arma::rowvec beta1, arma::mat P1) {
  // Dimensions
  int n = H.n_cols;
  int m = Q.n_cols;
  int t = H.n_slices;
  // Init matrices
  arma::mat a(t+1, m);
  arma::cube P(m, m, t+1);
  arma::mat nu(t, n);
  arma::cube F(n, n, t); 
  arma::cube Fi(n, n, t); 
  arma::cube K(m, n, t); 
  arma::cube L(m, m, t); 
  arma::mat r(t+1, m, arma::fill::zeros);
  arma::mat a_hat(t+1, m);
  // Prior for a1
  a.row(0) = beta1;
  P.slice(0) = P1;
  // Kalman filter
  for(int i = 0; i < t; i++) {
    nu.row(i) = y.row(i) - a.row(i) * Z.slice(i).t();
    F.slice(i) = Z.slice(i) * P.slice(i) * Z.slice(i).t() + H.slice(i);
    Fi.slice(i) = arma::inv(F.slice(i));
    K.slice(i) = Pi * P.slice(i) * Z.slice(i).t() * Fi.slice(i);
    L.slice(i) = Pi - K.slice(i) * Z.slice(i);
    a.row(i+1) = a.row(i) * Pi.t() + nu.row(i) * K.slice(i).t();
    P.slice(i+1) = Pi * P.slice(i) * L.slice(i).t() + Q.slice(i);
  }
  // Backward smoother
  for(int i = t-1; i >= 0; i--) {
    r.row(i) = nu.row(i) * Fi.slice(i).t() * Z.slice(i) + r.row(i+1) * L.slice(i);
  }
  //arma::rowvec r0 = nu.row(0) * Fi.slice(0).t() * Z.slice(0) + r.row(0) * L.slice(0);
  
  a_hat.row(0) = beta1 +  r.row(0) * P1.t();
  
  for(int i = 0; i < t; i++) {
    a_hat.row(i+1) = a_hat.row(i) * Pi.t()  + r.row(i+1) * Q.slice(i).t();
  }
  return a_hat;
}

// [[Rcpp::export]]
arma::mat rcppDrawBeta_DK(arma::mat y, arma::cube Z, arma::cube H, arma::cube Q, arma::mat Pi, arma::rowvec beta1, arma::mat P1) {
  // Dimensions
  int n = H.n_cols;
  int m = Q.n_cols;
  int t = H.n_slices;
  // Init matrices
  arma::mat e_plus(t, n);
  arma::mat u_plus(t, m);
  arma::mat a_plus(t+1, m);
  arma::mat y_plus(t, n);

  // Initial state for a+. Fix by Marek Jarocinski 2005
  a_plus.row(0) = rcppRmvnorm(1, arma::zeros(m), P1);
  
  for(int i = 0; i < t; i++) {
    e_plus.row(i) = rcppRmvnorm(1, arma::zeros(n), H.slice(i));
    u_plus.row(i) = rcppRmvnorm(1, arma::zeros(m), Q.slice(i));
    y_plus.row(i) = a_plus.row(i) * Z.slice(i).t() + e_plus.row(i);
    a_plus.row(i+1) = a_plus.row(i) * Pi.t() + u_plus.row(i);
  }
  
  arma::mat y_star = y - y_plus;
  
  // Kalman filter and disturbance smoother
  arma::mat a_star = rcppKalmanFilterSmooth(y_star, Z, H, Q, Pi, beta1, P1);
  
  return a_plus.rows(0, t-1) + a_star.rows(0, t-1);
}