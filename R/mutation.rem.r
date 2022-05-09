
#' Estimation of the random-effect model for clonality based on mutations.




mutation.rem <- function(mutmat, proba=FALSE, print.proba=FALSE, xigrid = seq(0.0005, 0.9995, by=0.001), init.para = c(0,1,0.5),
				conv.crit = 1e-5, niter=300){

  # Grid of likelihood values for each cases (rows) and each value of xi (columns)
  likmat <- grid.lik(xigrid, mutmat[,-1], mutmat[,1])
  # Maximum likelihood estimators
  out0 <- -100*abs(sum(log( init.para[3] * likmat[,-1] %*% xidens(init.para[1],init.para[2],xigrid)[-1]  * 0.001 + init.para[3] * likmat[,1] ) ))
  est <- clonEM(mutmat, init.para, xigrid, conv.crit, niter)
  out <- NULL
  out$mu <- est$para[1]
  out$sigma <- est$para[2]
  out$pi <- est$para[3]
  out$likelihood <- model.lik(est$para, likmat, out0, xigrid)
  out$convergence <- est$convergence
  out$n.iter <- est$n.iter
  out$likmat <- likmat
  
  if(proba==TRUE){
    pclon <- est$para[3] * rowSums(likmat[,-1] %*% xidens(est$para[1],est$para[2], xigrid)[-1] ) * 0.001
    pnclon <- (1-est$para[3]) * likmat[,1]
    out$pr.clonal <- pclon / (pclon+pnclon)
  }
  else{
    out$pr.clonal <- NULL
    print.proba <- FALSE
  }
  out$print.proba <- print.proba
  class(out) <- "mutation.rem"
  out
}

