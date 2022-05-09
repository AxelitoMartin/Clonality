oneseg <-
function(x, alpha, nperm, sbdry) {
# setup data and parameters for fortran code
  data.type <- "logratio"
  current.genomdat <- x - mean(x)
  current.tss <- sum(current.genomdat^2)
  n <- current.n <- length(x)
  sbn <- length(sbdry)
  
# there are used only by hybrid code but need to set it
  kmax <- 25
  nmin <- 200
  ngrid <- 100
  tol <- 1e-6

  if (n>200) {
    hybrid <- TRUE
    delta<-(kmax+1)/current.n
  } else {
    hybrid <- FALSE
    delta <- 0
  }

  # call the fortran segmentation cose
  zzz <- .Fortran("fndcpt",
                  n = as.integer(current.n),
                  x = as.double(current.genomdat),
                  tss = as.double(current.tss),
                  px = double(current.n),
                  sx = double(n),
                  nperm = as.integer(nperm),
                  cpval = as.double(alpha),
                  ncpt = integer(1),
                  icpt = integer(2),
                  ibin = as.logical(data.type=="binary"),
                  hybrid = as.logical(hybrid),
                  al0 = as.integer(2),
                  hk = as.integer(kmax),
                  delta = as.double(delta),
                  ngrid = as.integer(ngrid),
                  sbn = as.integer(sbn),
                  sbdry = as.integer(sbdry),
                  tol = as.double(tol),
                  PACKAGE = "DNAcopy")

# number of changes detected at alpha level (could be 0, 1 or 2)
  switch(1+zzz$ncpt, c(zzz$ncpt, 0, 0),
         c(zzz$ncpt, 1, zzz$icpt[1]),
         c(zzz$ncpt, c(1,0)+zzz$icpt))
}

