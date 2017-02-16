#'@importFrom Rcpp evalCpp
#'@useDynLib bayesVAR

### Classes

#' @export
setClass("bayesVAR_TVP")
#' @export
setClass("bayesFAVAR_TVP")
#' @export
setClass("bayesVAR")
#' @export
setClass("predictiveDensity_bayesVAR")

### Generics

#' Impulse response function for bayes VAR models
#' @export
impulse.response = function(x, ...) UseMethod("impulse.response", x)

#' Plot time-varying beta estimates for bayes VAR models
#' @export
plot.beta = function(x, ...) UseMethod("plot.beta", x)

#' Calculates predictive density for given bayes MCMC VAR model and number of periods
#' @param model bayes VAR model
#' @param T number of time steps ahead to predict
#' @export
predictive.density = function(x, ...) UseMethod("predictive.density", x)
