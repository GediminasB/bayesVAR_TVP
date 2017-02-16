#' Estimates FAVAR-TVP model using MCMC sampler.
#'
#' Estimates Time-Varying Parameters factor augmented VAR model using MCMC sampler.
#' \deqn{Y_t = F_t \Lambda + \nu_t; \nu_t ~ N(0, R)
#' \deqn{F_t = Z_t \beta_t + \epsilon_t; \epsilon_t ~ N(0, H)}
#' \deqn{\beta_t = \beta_{t-1} + u_t; u_t ~ N(0, Q)}
#' Prior parameters are estimated using OLS method on training sample.
#'
#' @param Y Matrix or timeseries object
#' @param l Number of time series to vectorize; all ts to vectorize must come before all other series
#' @param n.factors Number of factors tu use for first l series
#' @param p Integer for the lag order (default is p=1)
#' @param nburn Number of MCMC draws used to initialize the sampler (defaults to 10000)
#' @param nsim Number of MCMC draws excluding burn-in (defaults to 50000)
#' @param tau Length of training sample used for determining prior parameters via OLS
#' @param beta.algorithm Algorithm for drawing time-varying VAR parameters. Either 'DK' for Durbin and Coopman 2002, or 'CC' for Carter and Cohn 1994. Defaults to 'DK'
#' @param reject.explosive Should the explosive (not stable) draws of VAR parameters be rejected? Defaults to 'FALSE'. For rationale consult Cogley and Sargent 2005. May significantly increase computtion time!
#' @return \item{beta}{Draws of time-varying parameters (beta_t). 4D array [t x M x M*p+1 x draw]}
#'         \item{H}{Draws of observation equation error term covariance matrix H. 3D array [M x M x draw]}
#'         \item{Q}{Draws of state equation error term covariance matrix Q. 3D array [M x M x draw]}
#'         \item{L}{Draws of \Lambda. 3D array [n.factors x l x draw]}
#'         \item{F}{Draws of factors. 3D array }
#'         \item{R}{Draws of observation equation error term covariance matrix R.}
#' @references \itemize{
#'               \item \insertRef{KoopKorobilis2010}{bayesVAR}
#'               \item \insertRef{CC1994}{bayesVAR}
#'               \item \insertRef{DK2002}{bayesVAR}
#'               \item \insertRef{Bernanke2004}{bayesVAR}
#'             }
#' @export

#' @include bayesVAR_tvp.R
#' @export
bayesFAVAR_TVP = function(Y, l = ncol(Y), n.factors, p = 1, nburn = 10000, nsim = 50000, tau = 40, beta.algorithm = "DK", reject.explosive = FALSE) {
  # Prepare data
  Y_ = Y[-(1:tau),]
  # colnames(x.full) = c("const.", sapply(1:p, function(i) paste0(colnames(Y), "_L", i)))

  # Length vars
  N = nburn + nsim
  t = nrow(Y) - tau - p
  n = ncol(Y)
  n.f = n.factors+n-l
  n.vars = n.f*(n.f*p + 1)

  # Set the priors
  F.PCA = prcomp(Y[,1:l], center = FALSE, scale. = FALSE)
  F.PCA.x = F.PCA$x[, 1:n.factors]
  colnames(F.PCA.x) = paste0("PC", 1:n.factors)

  # Set the priors
  R.prior_scale = rep(0.01, l) # Bernanke et al. 2005
  R.prior_shape = 1*diag(l)  # Bernanke et al. 2005
  L.prior_M = diag(n.factors) # Bernanke et al. 2005
  L.prior_M.inv = solve(L.prior_M) # Bernanke et al. 2005

  # L.prior = t(F.PCA$rotation[,1:n.factors])
  F.prior = F.PCA.x[(tau+1):nrow(Y),]

  L.prior = solve(t(F.prior) %*% F.prior) %*% t(F.prior) %*% Y_[,1:l]
  rownames(L.prior) = colnames(F.PCA.x)

  # L.prior = L.prior/L.prior[1,1]
  # F.prior = F.prior*L.prior[1,1]

  # Rotate so that L[n.factors,n.factors] would be lower triangular
  # L.Rotation = solve(L.prior[,1:n.factors]) # diag(n.factors)

  U = chol(t(L.prior[,1:n.factors]) %*% L.prior[,1:n.factors])
  E = sqrt(solve(diag(diag(var(F.prior)))))
  L.Rotation = solve(t(U)) %*% t(L.prior[,1:n.factors])
  colnames(L.Rotation) = rownames(L.Rotation) = colnames(F.PCA.x)
  L.Rotation.inv = solve(L.Rotation)
  E = sqrt(solve(diag(diag(var(F.prior %*% L.Rotation.inv)))))

  L.prior = solve(E) %*% L.Rotation %*% L.prior
  L.prior[,1:n.factors][lower.tri(L.prior[,1:n.factors])] = 0

  F.prior = F.prior %*% L.Rotation.inv %*% E
  F.prior_P = 50 * diag(diag(var(F.prior)))

  # OLS prior. As in Koop Korobilis 2010. 4.1.1
  # We set aside part of data and use OLS solution for hyperparameters
  Y_OLS = cbind(F.PCA.x[1:tau,] %*% L.Rotation.inv, Y[1:tau,-(1:l)])
  y.train = Y_OLS[-(1:p),]
  x.train = cbind(`const.` = 1, lag(Y_OLS, 1:p)[-(1:p),])

  beta.OLS = t(solve(t(x.train) %*% x.train) %*% t(x.train) %*% y.train)
  resid.OLS = y.train - x.train %*% t(beta.OLS)
  sigma.OLS = 1/(tau - n.f*p - 1) * t(resid.OLS) %*% resid.OLS
  V.OLS = kronecker(solve(t(x.train) %*% x.train), sigma.OLS)

  # Residual variances of the corresponding p-lag univariate autoregressions
  # sigma.OLSi = apply(Y_OLS, 2, function(y) {
  #   X = as.matrix(lag.xts(y, 1:p), ncol = p)[-(1:p),]
  #   y_ = y[-(1:p)]
  #   beta.OLS = solve(t(X) %*% X) %*% t(X) %*% y_
  #   resid.OLS = y_ - X %*% beta.OLS
  #   var(resid.OLS)
  # })

  beta1.prior_mean = as.vector(matrix(beta.OLS, 1))
  beta1.prior_V = diag(diag(V.OLS)) # Ellis et al. 2014
  beta1.prior_V.inv = solve(beta1.prior_V)

  H.prior_nu = n.f + 1
  H.prior_S = diag(diag(sigma.OLS))

  Q.prior_nu = tau
  Q.prior_Q = 0.000001 * tau * V.OLS

  # Arrays in which we will keep draws
  beta.post = array(NA, c(t, n.vars, nsim))
  H.post = array(NA, c(n.f, n.f, nsim))
  Q.post = array(NA, c(n.vars, n.vars, nsim))
  F.post = array(NA, c(t+p, n.factors, nsim))
  L.post = array(NA, c(n.factors, l, nsim))
  R.post = array(NA, c(l, l, nsim))

  # Starting values
  H.draw = H.prior_S
  H.draw.inv = solve(H.draw)
  Q.draw = Q.prior_Q
  Q.draw.inv = solve(Q.draw)
  beta.draw = matrix(beta1.prior_mean, t, n.vars, byrow = TRUE)
  F.draw = F.prior
  L.draw = L.prior
  R.draw = diag(l)

  # Save some constants
  Q.post_nu = Q.prior_nu + t
  H.post_nu = H.prior_nu + t

  # Set up progress bar
  pb = progress::progress_bar$new(total = N, format = ":task | :current/:total [:bar]  :elapsed | @ :eta", clear = FALSE)
  pb$update(0, tokens = list(task = "Burn-in"))

  reject_n = 0; i = 2
  # Main loop
  while(i <= (N+1)) {
    F.draw_full = cbind(F.draw, Y[-(1:tau),-(1:l)])
    F.draw_full.x = cbind(1, lag(F.draw_full, 1:p))[-(1:p),]
    F.draw_full.y = F.draw_full[-(1:p),]
    Z = rcppExpandKronecker(F.draw_full.x, n.f)
    # Draw beta, contional on factors, H and Q
    beta1.post_V = chol2inv(chol(beta1.prior_V.inv + rcppZHZ(Z, H.draw.inv)))
    beta1.post_mean = beta1.post_V %*% (beta1.prior_V.inv %*% beta1.prior_mean + rcppZHy_TVP(Z, H.draw.inv, F.draw_full.y, beta.draw))
    beta.draw = .drawBeta(F.draw_full.y, Z, H.draw, Q.draw, beta1 = beta1.post_mean, P1 = beta1.post_V, algorithm = beta.algorithm)

    # Draw Q contional on factors and beta
    beta.draw_diff = beta.draw[-1,] - beta.draw[-t,]
    Q.post_Q.inv = chol2inv(chol(Q.prior_Q + t(beta.draw_diff) %*% beta.draw_diff))
    Q.draw.inv = rWishart(1, Q.post_nu, Q.post_Q.inv)[,,1]
    Q.draw = chol2inv(chol(Q.draw.inv))

    # Draw H conditional on data and beta
    H.post_S.inv = chol2inv(chol(H.prior_S + rcppSSEmat(F.draw_full.y, Z, beta.draw)))
    H.draw.inv = rWishart(1, H.post_nu, H.post_S.inv)[,,1]
    H.draw = chol2inv(chol(H.draw.inv))

    # Draw R conditional on factors and data
    Obs_OLS = chol2inv(chol(t(F.draw) %*% F.draw)) %*% t(F.draw) %*% Y_[,1:l]
    Obs_OLS.e = Y_[,1:l] - F.draw %*% Obs_OLS
    Obs_OLS.SSE = t(Obs_OLS.e) %*% Obs_OLS.e
    # Obs_SSE = Y_[,1:l] - F.draw %*% L.draw
    # R.bar = diag(t(Obs_SSE) %*% Obs_SSE) + R.prior_scale
    R.bar = Obs_OLS.SSE + diag(R.prior_scale) + t(Obs_OLS) %*% solve(L.prior_M.inv + solve(t(F.draw) %*% F.draw)) %*% Obs_OLS

    R.draw = diag(invgamma::rinvgamma(l, shape = t + p + diag(R.prior_shape), scale = diag(R.bar)))
    R.draw.inv = solve(R.draw)

    # Draw lambda on factors and data
    M.bar = t(F.draw) %*% F.draw + L.prior_M
    M.bar.inv = chol2inv(chol(M.bar))
    Gamma.bar = M.bar.inv %*% t(F.draw) %*% F.draw %*% Obs_OLS

    #Obs_OLS2 = chol2inv(chol(t(F.draw[,1:2]) %*% F.draw[,1:2])) %*% t(F.draw[,1:2]) %*% Y_[,2]
    #Gamma.bar2 = solve(M.bar[1:2,1:2]) %*% t(F.draw[,1:2]) %*% F.draw[,1:2] %*% Obs_OLS2

    L.draw = matrix(rcppRmvnorm(1, c(Gamma.bar[,-(1:(n.factors))]), kronecker(R.draw[-(1:(n.factors)),-(1:(n.factors))], M.bar.inv)), nrow = n.factors)
    # L.draw2 = matrix(rcppRmvnorm(1, Gamma.bar2, kronecker(R.draw[2,2], M.bar.inv[1:2,1:2])), nrow = n.factors-1)
    L.draw = cbind(L.prior[,1:(n.factors)], L.draw)

    # Draw factors conditional on data
    # ZHZ = nrow(Y_)*(L.draw %*% R.draw.inv %*% t(L.draw))
    # F1.post_V = chol2inv(chol(diag(n.factors) + ZHZ))
    # F1.post_mean = F1.post_V %*% (diag(n.factors) %*% rep(0, n.factors) + rcppZHy_FAVAR(t(L.draw), R.draw.inv, Y_[,1:l], F.draw))
    F.draw = .drawBeta(Y_[,1:l], t(L.draw), R.draw, H.draw[1:n.factors, 1:n.factors], beta1 = rep(0, n.factors), P1 = F.prior_P, algorithm = "CC")

    # Reject or save
    if(!reject.explosive || !.checkExplosive(beta.draw[t, ], n.f, p)) {
      if(i > nburn + 1) {
        j = i - nburn - 1
        beta.post[,,j] = beta.draw
        Q.post[,,j] = Q.draw
        H.post[,,j] = H.draw
        F.post[,,j] = F.draw
        L.post[,,j] = L.draw
        R.post[,,j] = R.draw
      }
      i = i + 1
    } else {
      reject_n = reject_n + 1
    }

    # Update the progress bar
    if((i-1) %% 100 == 0) pb$update((i-1)/N, tokens = list(task = ifelse(i-1 <= nburn, "Burn-in  ", "Sampling")))
  }
  # Add dimnames
  make.names_X.F = c("const.", sapply(1:p, function(i) paste0(colnames(Y_OLS), "_L", i)))
  make.names_betas = c(outer(colnames(Y_OLS), make.names_X.F, paste, sep = ":"))
  dim(beta.post) = c(t, n.f, n.f*p + 1, nsim)
  dimnames(beta.post)[[3]] = make.names_X.F
  dimnames(beta.post)[[2]] = colnames(Y_OLS)
  dimnames(beta.post)[[1]] = index(Y_)[-(1:p)]
  colnames(H.post) = rownames(H.post) = colnames(Y_OLS)
  colnames(Q.post) = rownames(Q.post) = make.names_betas
  colnames(R.post) = rownames(R.post) = colnames(L.draw) = colnames(Y)[1:l]
  rownames(L.draw) = colnames(Y_OLS)[1:n.factors]
  rownames(F.draw) = index(Y[-(1:tau),])
  colnames(F.draw) = colnames(Y_OLS)[1:n.factors]

  # Return
  structure(
    list(
      prior = list(
        beta0_mean = beta1.prior_mean,
        beta0_V = beta1.prior_V,
        H_nu = H.prior_nu,
        H_S = H.prior_S,
        Q_nu = Q.prior_nu,
        Q_Q = Q.prior_Q,
        F = F.prior,
        F_P = F.prior_P,
        L = L.prior,
        L.Rotation = L.Rotation,
        L_M = L.prior_M,
        R_scale = R.prior_scale,
        R_shape = R.prior_shape
      ),
      beta = beta.post,
      H = H.post,
      Q = Q.post,
      R = R.post,
      F = F.post,
      L = L.post,
      y = Y,
      factor.names = colnames(Y_OLS),
      var.names = colnames(Y_),
      t = t, n = n, n.vars = n.vars, p = p, nburn = nburn, nsim = nsim, l = l, n.factors = n.factors, n.f = n.f,
      rejection_rate = reject_n/(N+reject_n)
    ), class = "bayesFAVAR_TVP")
}

# Estimate of coefficients given loss function
#' @method coef bayesFAVAR_TVP
#' @export
coef.bayesFAVAR_TVP = function(model, loss.function = "quadratic") {
  if(loss.function == "quadratic") summary.func = mean
  if(loss.function == "0/1") summary.func = .mode.estimate
  if(loss.function == "absolute") summary.func = median
  beta.est = aperm(apply(model$beta, 1:3, summary.func), c(2, 3, 1))
  H.est = apply(model$H, 1:2, summary.func)
  Q.est = apply(model$Q, 1:2, summary.func)
  L.est = apply(model$L, 1:2, summary.func)
  F.est = apply(model$F, 1:2, summary.func)
  R.est = apply(model$R, 1:2, summary.func)

  list(beta.est = beta.est, H.est = H.est, Q.est = Q.est, L.est = L.est, F.est = F.est, R.est = R.est)
}

# Plot time varying parameters
#' @method plot.beta bayesFAVAR_TVP
#' @export
plot.beta.bayesFAVAR_TVP = plot.beta.bayesVAR_TVP

# Impulse response
#' @describeIn impulse.response Impulse response for bayes FAVAR-TVP
#' @method impulse.response bayesFAVAR_TVP
#' @export
impulse.response.bayesFAVAR_TVP = function(model, R = 20, t = model$t, orthogonal = TRUE, reorder = 1:model$n, plot = TRUE) {
  require(ggplot2)
  Phi = simplify2array(rcppIRF(model$beta[t,,-1,], model$H, R, orthogonal = orthogonal))
  dimnames(Phi) = list(impulse = model$factor.names, response = model$factor.names, t = (-model$p+1):R, path = NULL)

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

#' @describeIn predictive.density Predictive density for bayes TVP-VAR model.
#' @method predictive.density bayesFAVAR_TVP
#' @export
predictive.density.bayesFAVAR_TVP = function(model, T = 4, currentLevels = tail(model$y, model$p), n.sim = 10) {
  B = model$beta[model$t,,,]
  dim(B) = c(model$n.vars, model$nsim)
  Y.forecast = rcppPredictiveSim_FTVP(currentLevels, t(B), model$H, model$Q, model$L, model$R, T, k = n.sim)
  dimnames(Y.forecast) = list(t = (1-model$p):T, model$var.names, NULL)
  structure(Y.forecast, class = "predictiveDensity_bayesVAR")
}

