#' Generate Constraint Matrices
#' 
#' Generate a constraint matrix with a predefined structure
#' 
#' @param n a (possibly named) vector of sample sizes for each group
#' @param type character string defining the type of constraints; one of "monotone", "control","average","laverage","uaverage", or "caverage"
#' @param base column of the constraint matrix representing a control group (when type = "control")
#' 
#' @return a constraint matrix
#' 
#' @seealso \code{\link{orlm}}, \code{\link{constrSet}}
#' 
#' @keywords misc
#' 
#' @examples 
#' n <- c(10,20,30,40)
#' constrMat(n, type="monotone")
#' constrMat(n, type="control", base=2)
#' constrMat(n, type="average")
#' constrMat(n, type="laverage")
#' constrMat(n, type="uaverage")
#' constrMat(n, type="caverage", base=2)


constrMat <- function(n, type=c("monotone","control","average","laverage","uaverage","caverage"), base=1){
  d <- length(n)
  if (d < 2) stop("less than two groups")
  if (!is.numeric(n)) stop(sQuote("n"), " is not numeric")
  if (base < 1 || base > d) stop("base is not between 1 and ", d)
  if (!is.null(names(n))) varnames <- names(n) else varnames <- 1:d
  type <- match.arg(type)
  switch(type, monotone = {
    cm <- cbind(-diag(d-1),0) + cbind(0,diag(d-1))
  }, control = {
    cm <- diag(d)[-base,]
    cm[,base] <- -1
  }, average = {
    cm0 <- matrix(nrow=d-1, ncol=d)
    cm <- t(sapply(1:(d-1), function(i){
      cm0[i, 1:i] <- -n[1:i]/sum(n[1:i])
      cm0[i, (i+1):d] <- n[(i+1):d]/sum(n[(i+1):d])
      return(cm0[i,])
    }))
  }, laverage = {
    cm0 <- cbind(0,diag(d-1))
    cm <- t(sapply(1:(d-1), function(i){
      cm0[i, 1:i] <- -n[1:i]/sum(n[1:i])
      return(cm0[i,])
    }))
  }, uaverage = {
    cm0 <- cbind(-diag(d-1),0)
    cm <- t(sapply(1:(d-1), function(i){
      cm0[i, (i+1):d] <- n[(i+1):d]/sum(n[(i+1):d])
      return(cm0[i,])
    }))
  }, caverage = {
    cm0 <- matrix(0, nrow=d-1, ncol=d)
    cm0[,base] <- -1
    cm <- t(sapply(1:(d-1), function(i){
      id <- (1:d)[-base][1:i]
      cm0[i, id] <- n[id]/sum(n[id])
      return(cm0[i,])
    }))
  })
  colnames(cm) <- varnames
  return(cm)  
}
