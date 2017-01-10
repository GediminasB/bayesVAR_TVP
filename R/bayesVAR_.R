.makeZ = function(x, n) {
  t = nrow(x)
  n.vars = ncol(x) * n
  Z = array(NA, dim = c(n, n.vars, t))
  I = diag(n)
  colnames(I) = rownames(I) = colnames(x)[-1]
  for (i in 1:t) {
    Z[,,i] <- kronecker(x[i,], diag(n), make.dimnames = TRUE)
  }
  Z
}

.mode.estimate = function(x) {
  d = density(x)
  d$x[which.max(d$y)]
}

coef.bayesVAR = function(model, loss.function = "quadratic") {
  if(loss.function == "quadratic") summary.func = mean
  if(loss.function == "0/1”") summary.func = .mode.estimate
  if(loss.function == "absolute”") summary.func = median
  beta.est = apply(model$beta, 1:2, summary.func)
  H.inv.est = apply(model$H.inv, 1:2, summary.func)

  list(beta.est = beta.est, H.est = solve(H.inv.est))
}

bayesVAR = function(Y, nburn = 1000, nsim = 10000,
                    beta.prior_mean = rep(0, ncol(Y)*(ncol(Y) + 1)),
                    beta.prior_V = 0*diag(ncol(Y)*(ncol(Y) + 1)),
                    H.prior_nu = 2, H.prior_S = 0.0001*diag(ncol(Y))) {
  # require(mvtnorm)
  # require(progress)
  # require(Rcpp)
  # sourceCpp("rcppbayesVAR.cpp")

  # Inverse Prior matrices
  beta.prior_V.inv = beta.prior_V
  H.prior_S.inv = solve(H.prior_S)

  y = Y[-1,]
  x = na.omit(cbind(`const.` = 1, lag(Y)))

  N = nburn + nsim

  t = nrow(y)
  n = ncol(y)
  n.vars = n*ncol(x)

  Z = rcppExpandKronecker(x, n)

  beta.post = matrix(NA, nrow = N+1, ncol = n.vars)
  H.post.inv = array(NA, c(n, n, N+1))

  H.post.inv[,,1] = diag(n)
  beta.post[1,] = rep(0, n.vars)

  pb = progress::progress_bar$new(total = N, format = ":task | :current/:total [:bar]  :elapsed | @ :eta", clear = FALSE)
  pb$tick(0, tokens = list(task = "Burn-in  "))

  for(i in 2:(N+1)) {

    beta.post_V = solve(beta.prior_V.inv + rcppZHZ(Z, H.post.inv[,,i-1]))
    beta.post_mean = beta.post_V %*% (beta.prior_V.inv %*% beta.prior_mean + rcppZHy(Z, H.post.inv[,,i-1], y))
    beta.post[i,] = rmvnorm(1, beta.post_mean, beta.post_V)

    H.post_S.inv = solve(H.prior_S + rcppSSE(y, Z, beta.post[i,]))
    H.post_nu = H.prior_nu + t
    H.post.inv[,,i] = rWishart(1, H.post_nu, H.post_S.inv)[,,1]

    if((i-1) %% 100 == 0) pb$update((i-1)/N, tokens = list(task = ifelse(i-1 <= nburn, "Burn-in  ", "Sampling")))
  }
  b.out = beta.post[(nburn+2):(N+1),]
  dim(b.out) = c(nsim, n, n+1)
  b.out = aperm(b.out, c(2,3,1))
  colnames(b.out) = colnames(x)
  rownames(b.out) = colnames(y)

  H.out = H.post.inv[,,(nburn+2):(N+1)]
  colnames(H.out) = rownames(H.out) = colnames(y)

  structure(list(beta = b.out, H.inv = H.out), class = "bayesVAR")
}

# A = bayesVAR(data.selected[-1,])
#
# library(MTS)
# B = BVAR(as.ts(y), p = 1, C = 0.1*diag(5), V0 = diag(4))
