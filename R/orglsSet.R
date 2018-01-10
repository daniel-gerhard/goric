#' Set of generalised least-squares models
#' 
#' Fitting a specific set of generalisd least-squares models with order restrictions.
#' 
#' This function is just a wrapper for repeated calls of \code{\link{orgls}} with different constraint definitions. Predefined lists with constraint-sets can be constructed with function \code{\link{constrSet}}.
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param weights a \code{\link{varClasses}} object; more details are provided on the help pages in R package \code{nlme}
#' @param correlation a \code{\link{corClasses}} object; more details are provided on the help pages in R package \code{nlme}
#' @param set either a character string (see \code{\link{constrSet}}), or a list with slots for constr, rhs, and nec similarly defined as in \code{\link{orlm}}
#' @param direction direction of the order constraints
#' @param n a (possibly named) vector of sample sizes for each group
#' @param base column of the constraint matrix representing a control group
#' @param control a list of control arguments; see \code{\link{orlmcontrol}} for details.
#' 
#' @return a list with orgls objects
#' 
#' @seealso \code{\link{orgls}}, \code{\link{constrSet}}, \code{\link{goric}}
#' 
#' @keywords models


orglsSet <- function(formula, data, weights=NULL, correlation=NULL, set, direction="increase", n=NULL, base=1, control=orlmcontrol()){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf, contrasts)
  d <- ncol(x)
  if (is.null(n)){
    n <- rep(1,d)
    names(n) <- colnames(x)
  }
  if (length(n) != d) stop("n has not the same length as there are columns in the design matrix!")
  if (is.character(set)) set <- constrSet(n, set=set, direction=direction, base=base)
  out <- lapply(set, function(s){
    orgls(formula, data, constr=s$constr, rhs=s$rhs, nec=s$nec, control=control, weights=weights, correlation=correlation)
  })
  class(out) <- c("list","orglslist")
  return(out)
}