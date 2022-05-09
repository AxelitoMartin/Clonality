#____________________       Auxiliary functions


clonEM <- function(mutmat, init.para, xigrid, conv.crit, niter) {
 
  probamut <- mutmat[, 1]
  mutns <- mutmat[, -1]
  # number of cases
  ncase <- ncol(mutns)
  # number of grid points
  ngrid <- length(xigrid)
  # initialize parameters (mu, sigma, pi)
  para <- init.para
  # probability of matches for given clonality signal
  pYj0 <- grid.lik(0, mutns, probamut)
  pYj <- grid.lik(xigrid, mutns, probamut)
  # clonality signal density for xigrid
  g <- xidens(para[1], para[2], xigrid)
  gj <- rep(0, ngrid)
  vxij <- exij <- vpij <- rep(NA_real_, ncase)
  zz <- log(-log(1-xigrid))
  for(j in 1:ncase) {
    # product density
    gj <- pYj[j,]*g
    # integral
    igj <- sum(gj)/ngrid
    # conditional estimate of case being clonal
    vpij[j] <- para[3]*igj/((1-para[3])*pYj0[j] + para[3]*igj)
    # expected value of zz 
   exij[j] <- sum(zz*gj)/sum(gj)
    # variance of zz 
    vxij[j] <- sum((zz^2)*gj)/sum(gj) - exij[j]^2
  }
  # estimate the parameters from maximizing the conditional expectation
  para[1] <- sum(vpij*exij)/sum(vpij)
  para[2] <- max(0.1, sqrt(sum(vpij*(exij - para[1])^2)/sum(vpij)))
  para[3] <- sum(vpij)/ncase
  # iterate through the parameters
  ii <- 1
  conv <- 0
  while(conv == 0 & ii <= niter) {
    prev <- para
    g <- xidens(para[1], para[2], xigrid)
    for(j in 1:ncase) { 
      gj <- pYj[j,]*g
      igj <- sum(gj)/ngrid
      vpij[j] <- para[3]*igj/((1-para[3])*pYj0[j] + para[3]*igj)
      exij[j] <- sum(zz*gj)/sum(gj)
    }
    para[1] <- sum(vpij*exij)/sum(vpij)
    # variance is expectation of variance + variance of expectation
    vxi <- sum(vpij*vxij)/sum(vpij) + sum(vpij*(exij - para[1])^2)/sum(vpij)
    # lower bound sigma by 0.1
    para[2] <- max(0.1, sqrt(vxi))
    para[3] <- sum(vpij)/ncase
     #convergence reached
    if( max(abs(para - prev)) < conv.crit ) conv <- 1
    ii <- ii + 1
 }
  return( list(para = para, n.iter = (ii-1), convergence = conv) )
}

