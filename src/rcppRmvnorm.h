#ifndef PKG_BayesVAR_H
#define PKG_BayesVAR_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
arma::mat rcppRmvnorm(int n, arma::colvec mu, arma::mat sigma);

#endif