#____________________       Auxiliary functions


#' Auxiliary function computing the density of xi
#'
#' @description Density function for the random variable xi, using a lognormal density for phi=-log(1-xi)
#' @param pmu Mean parameter of the distribution.
#' @param psig Variance parameter of the distribution.
#' @param xigrid Grid of the values of xi, corresponding to its domain of definition. 
#' @return Returns the density value for the given values of xi.
#' @export

xidens <- function(pmu, psig, xigrid) { sapply(xigrid, function(x) dlnorm(-log(1-x), pmu, psig)/(1-x)) }

