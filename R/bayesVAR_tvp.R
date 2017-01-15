#'@importFrom Rcpp evalCpp
#'@useDynLib bayesVAR

#' @export
setClass("bayesVAR_TVP")

#' @export
setClass("bayesVAR")

#' @export
setClass("predictiveDensity_bayesVAR")

# Impulse response generic
#' Impulse response function for bayes VAR models
#' @export
impulse.response = function(x, ...) UseMethod("impulse.response", x)

# Plot beta generic
#' Plot time-varying beta estimates for bayes VAR models
#' @export
plot.beta = function(x, ...) UseMethod("plot.beta", x)

# Predictive density generic
#' Calculates predictive density for given bayes MCMC VAR model and number of periods
#' @param model bayes VAR model
#' @param T number of time steps ahead to predict
#' @export
predictive.density = function(x, ...) UseMethod("predictive.density", x)

# Check if VAR parameters are explosive
# TODO: speed comparison with Arma implementation
.checkExplosive = function(beta, n, o) {
  F.upper = matrix(beta, nrow = n, ncol = n*p+1)[,-1]
  F.lower = cbind(diag((p-1)*n), matrix(0, nrow = (p-1)*n, ncol = n))
  F_ = rbind(F.upper, F.lower)
  any(abs(eigen(F_)$values) > 1)
}

# Conditional draw of beta. C++ function wrapper
.drawBeta = function(y, Z, H, Q, Pi = diag(ncol(Q)), beta1, P1 = diag(ncol(Q)), algorithm = "DK") {
  if(class(H) == "matrix") H = array(H, dim = c(nrow(H), ncol(H), nrow(y)))
  if(class(Q) == "matrix") Q = array(Q, dim = c(nrow(Q), ncol(Q), nrow(y)))
  if(algorithm == "DK") rcppDrawBeta = rcppDrawBeta_DK
  if(algorithm == "CC") rcppDrawBeta = rcppDrawBeta_CC
  rcppDrawBeta(y, Z, H, Q, Pi, beta1, P1)
}

#' Estimates TVP-VAR model using MCMC sampler.
#'
#' Estimates Time-Varying Parameters VAR model using MCMC sampler.
#' \deqn{Y_t = Z_t \beta_t + \epsilon_t; \epsilon_t ~ N(0, H)}
#' \deqn{\beta_t = \beta_{t-1} + u_t; u_t ~ N(0, Q)}
#' Prior parameters are estimated using OLS method on training sample.
#'
#' @param Y Matrix or timeseries object
#' @param p Integer for the lag order (default is p=1)
#' @param nburn Number of MCMC draws used to initialize the sampler (defaults to 10000)
#' @param nsim Number of MCMC draws excluding burn-in (defaults to 50000)
#' @param tau Length of training sample used for determining prior parameters via OLS
#' @param beta.algorithm Algorithm for drawing time-varying VAR parameters. Either 'DK' for Durbin and Coopman 2002, or 'CC' for Carter and Cohn 1994. Defaults to 'DK'
#' @param reject.explosive Should the explosive (not stable) draws of VAR parameters be rejected? Defaults to 'FALSE'. For rationale consult Cogley and Sargent 2005. May significantly increase computtion time!
#' @return \item{beta}{Draws of time-varying parameters (beta_t). 4D array [t x M x M*p+1 x draw]}
#'         \item{H}{Draws of observation equation error term covariance matrix H. 3D array [M x M x draw]}
#'         \item{Q}{Draws of state equation error term covariance matrix Q. 3D array [M x M x draw]}
#' @references \itemize{
#'               \item \insertRef{KoopKorobilis2010}{bayesVAR}
#'               \item \insertRef{CC1994}{bayesVAR}
#'               \item \insertRef{DK2002}{bayesVAR}
#'               \item \insertRef{CogleySargent2005}{bayesVAR}
#'             }
#' @export
bayesVAR_TVP = function(Y, p = 1, nburn = 10000, nsim = 50000, tau = 40, beta.algorithm = "DK", reject.explosive = FALSE) {
  # Prepare data
  y.full = Y[-(1:p),]
  # x.full = na.omit(cbind(`const.` = 1, Y[-nrow(Y),]))
  x.full = cbind(1, lag(Y, 1:p))[-(1:p),]
  colnames(x.full) = c("const.", sapply(1:p, function(i) paste0(colnames(Y), "_L", i)))

  y.train = y.full[1:tau,]
  x.train = x.full[1:tau,]
  y = y.full[(tau+1):nrow(y.full),]
  x = x.full[(tau+1):nrow(x.full),]

  # Length vars
  N = nburn + nsim
  t = nrow(y)
  n = ncol(y)
  n.vars = n*ncol(x)

  # OLS prior. As in Koop Korobilis 2010. 4.1.1
  # We set aside part of data and use OLS solution for hyperparameters
  beta.OLS = t(solve(t(x.train) %*% x.train) %*% t(x.train) %*% y.train)
  resid.OLS = y.train - x.train %*% t(beta.OLS)
  sigma.OLS = 1/(tau - n*p - 1) * t(resid.OLS) %*% resid.OLS
  V.OLS = kronecker(solve(t(x.train) %*% x.train), sigma.OLS)
  # Set the priors
  beta0.prior_mean = matrix(beta.OLS, 1)
  beta0.prior_V = 4 * V.OLS
  H.prior_nu = n + 1
  H.prior_S = diag(n)
  Q.prior_nu = tau
  Q.prior_Q = 0.0001 * tau * V.OLS

  # Expand x into 3d array
  Z = rcppExpandKronecker(x, n)
  # Arrays in which we will keep draws
  beta.post = array(NA, c(t, n.vars, N+1))
  H.post = array(NA, c(n, n, N+1))
  Q.post = array(NA, c(n.vars, n.vars, N+1))
  # Starting values
  H.draw = H.prior_S
  H.post.inv = solve(H.post[,,1])
  Q.draw = Q.prior_Q
  Q.post.inv = solve(Q.post[,,1])
  beta.post[,,1] = matrix(0, t, n.vars)
  # Set up progress bar
  pb = progress::progress_bar$new(total = N, format = ":task | :current/:total [:bar]  :elapsed | @ :eta", clear = FALSE)
  pb$update(0, tokens = list(task = "Burn-in"))

  reject_n = 0; i = 2
  # Main loop
  while(i <= (N+1)) {
    # Draw beta, contional on data, H and Q
    beta.draw = .drawBeta(y, Z, H.draw, Q.draw, beta1 = beta0.prior_mean, P1 = beta0.prior_V, algorithm = beta.algorithm)
    # Draw Q contional on data and beta
    beta.post_diff = beta.draw[-1,] - beta.draw[-t,]
    Q.post_Q.inv = chol2inv(chol(Q.prior_Q + t(beta.post_diff) %*% beta.post_diff))
    Q.post_nu = Q.prior_nu + t
    Q.post.inv = rWishart(1, Q.post_nu, Q.post_Q.inv)[,,1]
    Q.draw = chol2inv(chol(Q.post.inv))
    # Draw H conditional on data and beta
    H.post_S.inv = chol2inv(chol(H.prior_S + rcppSSEmat(y, Z, beta.draw)))
    H.post_nu = H.prior_nu + t
    H.post.inv = rWishart(1, H.post_nu, H.post_S.inv)[,,1]
    H.draw = chol2inv(chol(H.post.inv))

    if(!reject.explosive || !.checkExplosive(beta.draw[t, ], n, p)) {
      beta.post[,,i] = beta.draw
      Q.post[,,i] = Q.draw
      H.post[,,i] = H.draw
      i = i + 1
    } else {
      reject_n = reject_n + 1
    }

    # Update the progress bar
    if((i-1) %% 100 == 0) pb$update((i-1)/N, tokens = list(task = ifelse(i-1 <= nburn, "Burn-in  ", "Sampling")))
  }
  # Subset the burnin and add dimnames
  b.out = beta.post[,, (nburn+2):(N+1)]
  dim(b.out) = c(t, n, n*p + 1, nsim)
  dimnames(b.out)[[3]] = colnames(x)
  dimnames(b.out)[[2]] = colnames(y)
  dimnames(b.out)[[1]] = index(y)
  H.out = H.post[,,(nburn+2):(N+1)]
  colnames(H.out) = rownames(H.out) = colnames(y)
  Q.out = Q.post[,,(nburn+2):(N+1)]
  colnames(Q.out) = rownames(Q.out) = colnames(Z)
  # Return
  structure(list(beta = b.out, H = H.out, Q = Q.out, y = y,
                 var.names = colnames(y), t = t, n = n, n.vars = n.vars, p = p, nburn = nburn, nsim = nsim, rejection_rate = reject_n/(N+reject_n)), class = "bayesVAR_TVP")
}

#' @describeIn predictive.density Predictive density for bayes TVP-VAR model.
#' @method predictive.density bayesVAR_TVP
#' @export
predictive.density.bayesVAR_TVP = function(model, T = 10) {
  B = model$beta[model$t,,,]
  dim(B) = c(model$n.vars, model$nsim)
  Y.forecast = rcppPredictiveSim_TVP(tail(model$y,1), t(B), model$H, model$Q, 10)
  dimnames(Y.forecast) = list(t = 0:T, model$var.names, NULL)
  structure(Y.forecast, class = "predictiveDensity_bayesVAR")
}

#' @method plot predictiveDensity_bayesVAR
#' @export
plot.predictiveDensity_bayesVAR = function(Y.forecast, w = NULL) {
  require(ggplot2)
  Y.fc_quantile = apply(Y.forecast, 1:2, Hmisc::wtd.quantile, probs = c(0.05, 0.16, 0.5, 0.84, 0.95), weights = w, normwt = TRUE)
  dimnames(Y.fc_quantile)[[1]] = paste0("Q_", c(0.05, 0.16, 0.5, 0.84, 0.95)*100)
  Y.fc_quantile.melt = reshape2::melt(Y.fc_quantile, varnames = c("quantile", "t", "ts"))
  Y.fc_quantile.dcast = reshape2::dcast(Y.fc_quantile.melt, ts + t ~ quantile)
  ggplot(Y.fc_quantile.dcast, aes(x = t)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_ribbon(aes(ymax = Q_95, ymin = Q_5, fill = "90%"), alpha = .25) +
    geom_ribbon(aes(ymax = Q_84, ymin = Q_16, fill = "68%"), alpha = .25) +
    geom_line(aes(color = "median", y = Q_50)) +
    facet_wrap(~ts) +
    scale_colour_manual("", values = "black") +
    scale_fill_manual("", values = c("grey12", "black")) +
    theme_bw() + theme(legend.direction = "horizontal", legend.position = "top")
}

#' Minimum entropy weights
#'
#' Finds minimum entropy weights using Kullback-Leibler information criterion.
#' New prediction density conditions are expressed by function \eqn{g} and vector \eqn{g.rhs} such that
#' \deqn{\mathbb{E}g(pd) = g.rhs}{Eg(pd) = g.rhs}
#'
#' @param pd Simulations of predictive density [T+1 x n x N.sim]
#' @param g Transformation for pd used to express to new information. Must return vector of length p
#' @param g.rhs Vector of length p
#' @return \item{pi.star}{Adjusted probabilities}
#'         \item{KLIC}{Kullback-Leibler information criterion}
#'         \item{gamma}{Vector of Lagrange multipliers}
#'         \item{optim.code}{An integer indicating why the optimization process terminated. See ?nlm for info}
#' @references \itemize{
#'               \item \insertRef{Robertson2005}{bayesVAR}
#'             }
#'
#' @export
MinimumEntropy_weights = function(pd, g, g.rhs) {
  N = dim(pd)[3]
  Y.transf = apply(pd, 3, g)
  Y.transf.dist = Y.transf - g.rhs
  f.min = function(gamma) mean(exp(t(gamma) %*%  Y.transf.dist))
  f.min.start = -rowSums(Y.transf.dist)/rowSums(Y.transf.dist^2) # using Maclaurin series expansion of exp at x = 0 (e^(c*x) = 1 + c*x)
  f.min.solve = stats::nlm(f.min, p = f.min.start)
  gamma = f.min.solve$estimate
  pi.star_ = exp(t(gamma) %*% Y.transf)[1,]
  pi.star = pi.star_/sum(pi.star_)
  KLIC = sum(pi.star * (log(pi.star) + log(N)))

  list(pi.star = pi.star, KLIC = KLIC, gamma = gamma, optim.code = f.min.solve$code)
}


# Estimate of coefficients given loss function
#' @method coef bayesVAR_TVP
#' @export
coef.bayesVAR_TVP = function(model, loss.function = "quadratic") {
  if(loss.function == "quadratic") summary.func = mean
  if(loss.function == "0/1”") summary.func = .mode.estimate
  if(loss.function == "absolute”") summary.func = median
  beta.est = aperm(apply(model$beta, 1:3, summary.func), c(2, 3, 1))
  H.est = apply(model$H, 1:2, summary.func)
  Q.est = apply(model$Q, 1:2, summary.func)

  list(beta.est = beta.est, H.est = H.est, Q.est = Q.est)
}

# Plot time varying parameters
#' @method plot.beta bayesVAR_TVP
#' @export
plot.beta.bayesVAR_TVP = function(model) {
  require(ggplot2)
  A.coef = coef(model)
  A.coef.melt = reshape2::melt(A.coef$beta.est, varnames = c("row", "col", "t"))
  A.plot = ggplot(A.coef.melt, aes(x = t, y = value)) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    geom_line() + facet_grid(col ~ row, scales = "free_y") + theme_bw()
  print(A.plot)
}

# Impulse response
#' @describeIn impulse.response Impulse response for bayes VAR-TVP
#' @method impulse.response bayesVAR_TVP
#' @export
impulse.response.bayesVAR_TVP = function(model, R = 20, t = model$t, orthogonal = TRUE, reorder = 1:model$n, plot = TRUE) {
  require(ggplot2)
  Phi = simplify2array(rcppIRF(model$beta[t,,-1,], model$H, R, orthogonal = orthogonal))
  dimnames(Phi) = list(impulse = model$var.names, response = model$var.names, t = (-model$p+1):R, path = NULL)

  Phi.q = apply(Phi, 1:3, quantile, c(0.05, 0.16, 0.5, 0.84, 0.95))
  Phi.q.melt = reshape2::melt(Phi.q[,,,-1], varnames = c("quantile", "impulse", "response", "t"))
  Phi.q.dcast = reshape2::dcast(Phi.q.melt, impulse + response + t ~ quantile)

  if(plot)
    print(ggplot(Phi.q.dcast, aes(x = t)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_ribbon(aes(ymax = `95%`, ymin = `5%`, fill = "90%"), alpha = .25) +
      geom_ribbon(aes(ymax = `84%`, ymin = `16%`, fill = "68%"), alpha = .25) +
      geom_line(aes(y = `50%`, color = "Median")) +
      facet_grid(impulse~response) +
      xlab("Response") + ylab("Impulse") +
      scale_colour_manual("", values = "black") +
      scale_fill_manual("", values = c("grey12", "black")) +
      theme_bw() + theme(legend.direction = "horizontal", legend.position = "top") +
      coord_cartesian(ylim = c(-1, 1)))

  # Return
  Phi
}
