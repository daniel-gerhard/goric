#' Control arguments for the orlm function.
#' 
#' A list with control arguments controlling the orlm function
#' 
#' @param maxiter maximum number of iterations
#' @param absval tolerance criterion for convergence
#' 
#' @return a list with control arguments
#' 
#' @seealso \code{\link{orlm}}
#' 
#' @keywords models

orlmcontrol <-
function(maxiter=10000, absval=0.0001){
  if (maxiter < 1) stop("max No of iterations < 1")
  if (absval < .Machine$double.eps){
    warning(paste("absval <", .Machine$double.eps,"\nabsval is set to", .Machine$double.eps))
    absval <- .Machine$double.eps
  }
  return(list(maxiter=maxiter, absval=absval))
}

