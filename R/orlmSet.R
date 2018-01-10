#' Set of multivariate regression models
#' 
#' Fitting a specific set of multivariate regression models with order restrictions.
#' 
#' This function is just a wrapper for repeated calls of \code{\link{orlm}} with different constraint definitions. Predefined lists with constraint-sets can be constructed with function \code{\link{constrSet}}.
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param set either a character string (see \code{\link{constrSet}}), or a list with slots for constr, rhs, and nec similarly defined as in \code{\link{orlm}}
#' @param direction direction of the order constraints
#' @param n a (possibly named) vector of sample sizes for each group
#' @param base column of the constraint matrix representing a control group
#' @param control a list of control arguments; see \code{\link{orlmcontrol}} for details.
#' 
#' @return a list with orlm objects
#' 
#' @seealso \code{\link{orlm}}, \code{\link{constrSet}}, \code{\link{goric}}
#' 
#' @keywords models
#' 
#' @examples 
#' ########################
#' ## Artificial example ##
#' ########################
#' 
#' n <- 10
#' m <- c(1,2,4,5,2,1)
#' nm <- length(m)
#' dat <- data.frame(grp=as.factor(rep(1:nm, each=n)),
#'                   y=rnorm(n*nm, rep(m, each=n), 1))
#' 
#' (cs <- constrSet(table(dat$grp), set="sequence"))
#' (oss <- orlmSet(y ~ grp-1, data=dat, set=cs))
#' 
#' # the same as:
#' oss <- orlmSet(y ~ grp-1, data=dat, set="sequence")

orlmSet <- function(formula, data, set, direction="increase", n=NULL, base=1, control=orlmcontrol()){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
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
    orlm(formula, data, constr=s$constr, rhs=s$rhs, nec=s$nec, control=control)
  })
  class(out) <- c("list","orlmlist")
  return(out)
}


