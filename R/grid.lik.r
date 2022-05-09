#____________________       Auxiliary functions


#' Auxiliary function: Grid of conditional probabilities
#'
#' This auxiliary function generates the grid of likelihood values for each tumor pair (rows) and each value of xi (columns): P(observed mutations | xi)
#' @param xigrid Grid of the values of xi, corresponding to its domain of definition. 
#' @param mutns Matrix of the mutations observed, with all mutations in rows and the cases (tumor pairs) in columns. 
#' The data are coded as 0=mutation not observed, 1=shared mutation (observed in both tumors), 2=private mutation 
#' (observed in one tumor only). 
#' @param probamut Vector of the probabilities of occurence for each mutation.
#' @return Return the matrix of the likelihood values for each tumor pair (rows) and each value of xi (columns). This matrix is called 
#' by the auxiliary function grid.lik, returned as a parameter by the function clonal.est, and used as a parameter by the function clonal.proba.
#' @export


grid.lik <- function(xigrid, mutns, probamut) {
  # find the index of shared and private mutations for each case
  mutidx <- apply(mutns, 2, function(muts) {
    ishared <- which(muts==1)
    iprivate <- which(muts==2)
    list(ishared, iprivate)
  })
  sapply(xigrid, function(xi, probamut, mutidx) {
    # p is the conditional prob of shared mutation
    p <- {xi + {1-xi}*probamut}/{xi + {1-xi}*{2-probamut}}
    # helper vectors
    logp <- log(p)
    log1p <- log(1-p)
    
    unlist(lapply(mutidx, function(ii, logp, log1p) {
      exp(sum(logp[ii[[1]]]) + sum(log1p[ii[[2]]]))
    }, logp, log1p))
  }, probamut, mutidx)
}
