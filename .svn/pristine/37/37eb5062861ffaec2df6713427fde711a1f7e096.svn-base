
#' Print for the mutation.proba function
#'
#' @description Print a summary of results for the probabilities of clonality estimated by the mutation.proba function.  
#' @usage ## S3 method for class 'mutation.proba'
#' print(x, ...)
#' @param x a mutation.proba object
#' @param ... Other unused arguments.
#' @return Print results for the individual probabilities of clonality.
#' @seealso mutation.proba
#' @export
#' @method print mutation.proba
#' @S3method print mutation.proba

print.mutation.proba <- function(x, ...) 
{
  
  if(class(x)!="mutation.proba"){
    stop("The object x must be a class mutation.proba.")
  }else{
    if(!is.null(names(x))){
      outp <- cbind.data.frame(names(x), round(as.vector(x),3))
      colnames(outp) <- c("Tumor pairs", "Probability of being clonal")
    }
    else{
      outp <-  round(as.vector(x),3)
    }
    print(outp, row.names = F)
  }
}
