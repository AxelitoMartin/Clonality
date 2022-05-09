# SNVtest2 is similar to SNVtest except that
# pfreq is a 2-column matrix giving site specific mutation probabilities
SNVtest2 <- function (tumor1, tumor2, pfreq, nrep = 999) {
  # number of potentially mutated loci
  n <- nrow(pfreq)
  # consistency checks
  if (n != length(tumor1) | n != length(tumor2)) 
    stop("length of tumor1 and tumor2 should be equal nrow of pfreq")
  if (any(!(tumor1 %in% c(0, 1))) | any(!(tumor2 %in% c(0, 1)))) 
    stop("tumor1 and tumor2 should be a vectors with only values 0 and 1")
  if (any(pfreq > 1 | pfreq < 0)) 
    stop("pfreq should take numerical values between 0 and 1")
  
  # site sepcific probabilities
  pprob <- pfreq[,1] # site 1
  qprob <- pfreq[,2] # site 2
  # mutated loci
  ii <- which(tumor1 + tumor2 > 0)
  
  # generate mutations under independent tumor model
  cn <- length(ii)
  mprob1 <- pprob[ii]*qprob[ii]/(1-(1-pprob[ii])*(1-qprob[ii]))
  mprob2 <- mprob1 + pprob[ii]*(1-qprob[ii])/(1-(1-pprob[ii])*(1-qprob[ii]))
  # generate uniform numbers
  uu <- runif(nrep*cn)
  # initialize all as type 3
  mtype <- rep(3, nrep*cn)
  # if uu less than mprob1 subtract 1
  mtype <- mtype - 1*(uu < rep(mprob1, nrep))
  # if uu less than mprob2 subtract 1
  mtype <- mtype - 1*(uu < rep(mprob2, nrep))
  # so if (0 < uu < mprob1) type is 1 and (mprob1 < uu < mprob2) it is 2
  
  # stack observed and simulated mutation data together
  omtype <- rep(3, cn)
  omtype[tumor1[ii]==1 & tumor2[ii]==1] <- 1
  omtype[tumor1[ii]==1 & tumor2[ii]==0] <- 2
  mutns <- cbind(omtype, matrix(mtype, ncol=nrep))
  
  # likelihood function
  Lik <- function(ksi, x, p) {
    a <- which(x == 1)
    b <- which(x == 0)
    e <- c(a, b)
    sum(log((ksi * p + (1 - ksi) * p^2)[a])) + sum(log(2 * 
                                                         (1 - ksi) * (p * (1 - p))[b])) - sum(log((ksi * p + 
                                                                                                     (1 - ksi) * p^2 + 2 * (1 - ksi) * p * (1 - p))[e]))
  }
  # optimization function for clonality signal
  fLR <- function(x, p) {
    o <- optim(0.5, function(ksi) {
      Lik(ksi, x, p)
    }, control = list(fnscale = -1), method = "L-BFGS-B", 
    #    lower = 1e-05, upper = 0.999999)
    lower = 1e-05, upper = 0.999999)
    o$par
  }
  
  # estimated clonality signal for the two sites
  xihat1 <- apply(1*(mutns==1), 2, fLR, pprob[ii])
  xihat2 <- apply(1*(mutns==1), 2, fLR, qprob[ii])
  
  # sites 1 & 2 likelihood-ratios
  s12 <- apply(rbind(xihat1, xihat2, mutns), 2, function(x, p1, q1) {
    xi1 <- x[1]
    xi2 <- x[2]
    muts <- x[-(1:2)]
    ii1 <- which(muts==1)
    ii2 <- which(muts==2)
    ii3 <- which(muts==3)
    # site 1 likelihood ratio
    llr1 <- sum(log({xi1 + (1-xi1)*p1}/{p1*q1})[ii1])
    llr1 <- llr1 + sum(log({(1-xi1)*(1-p1)}/{p1*(1-q1)})[ii2])
    llr1 <- llr1 + sum(log({(1-xi1)*(1-p1)}/{(1-p1)*q1})[ii3])
    llr1 <- llr1 - sum(log({xi1 + (1-xi1)*(2-p1)}/{p1 + q1 - p1*q1}))
    # site 2 likelihood ratio
    llr2 <- sum(log({xi2 + (1-xi2)*q1}/{p1*q1})[ii1])
    llr2 <- llr2 + sum(log({(1-xi2)*(1-q1)}/{p1*(1-q1)})[ii2])
    llr2 <- llr2 + sum(log({(1-xi2)*(1-q1)}/{(1-p1)*q1})[ii3])
    llr2 <- llr2 - sum(log({xi2 + (1-xi2)*(2-q1)}/{p1 + q1 - p1*q1}))
    c(llr1, llr2)
  }, pprob[ii], qprob[ii])
  
  lrmax <- apply(s12, 2, max)
  pval <- pmax(sum(lrmax >= lrmax[1])/(nrep+1),1/(nrep+1))
  out <- c(sum(omtype==1), sum(omtype==2), sum(omtype==3), xihat1[1], xihat2[1],  pval)
  names(out) <- c("n.match", "n.site1only", "n.site2only", "xi.site1", "xi.site2", "p.value")
  out
}