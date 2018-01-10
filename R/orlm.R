#' Fitting multivariate regression models with order restrictions
#' 
#' This is a modification of the \code{lm} function, fitting (multivariate) linear models with order constraints on the model coefficients.
#' 
#' The contraints in the hypothesis of interest are defined by \eqn{Constr}, \eqn{rhs}, and \eqn{nec}. The first \eqn{nec} constraints are the equality contraints: \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}; and the remaing ones are the inequality contraints: \eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}.
#' Two requirements should be met:
#' \enumerate{
#'   \item The first \eqn{nec} constraints must be the equality contraints (i.e., \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}) and the remaining ones the inequality contraints (i.e., \eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}).
#'   \item When \eqn{rhs} is not zero, \eqn{Constr} should be of full rank (after discarding redundant restrictions).
#' } 
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param constr matrix with constraints; with rows as constraint definition, columns should be in line with the parameters of the model
#' @param rhs vector of right hand side elements; \eqn{Constr \; \theta \geq rhs}; number should equal the number of rows of the constr matrix
#' @param nec number of equality constraints; a numeric value treating the first nec constr rows as equality constraints, or a logical vector with \code{TRUE} for equality- and \code{FALSE} for inequality constraints.
#' @param control a list of control arguments; see \code{\link{orlmcontrol}} for details.
#' 
#' @return an object of class orlm
#' 
#' @references 
#' \itemize{
#' \item Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, \bold{98}, 495--501.
#' \item Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the Order-Restricted Information Criterion for Multivariate Normal Linear Models. \emph{Journal of Statistical Planning and Inference}, \bold{142}, 2454-2463. doi:10.1016//j.jspi.2012.03.007.
#' \item Kuiper R.M. and Hoijtink H. (submitted). A Fortran 90 Program for the Generalization of the Order-Restricted Information Criterion. Journal of Statictical Software.
#' }
#' 
#' @seealso \code{\link{solve.QP}}, \code{\link{goric}}
#' 
#' @keywords models
#'
#' @examples 
#' ########################
#' ## Artificial example ##
#' ########################
#' n <- 10
#' m <- c(1,2,1,5)
#' nm <- length(m)
#' dat <- data.frame(grp=as.factor(rep(1:nm, each=n)),
#'                   y=rnorm(n*nm, rep(m, each=n), 1))
#'
#' # unrestricted linear model
#' cm1 <- matrix(0, nrow=1, ncol=4)
#' fm1 <- orlm(y ~ grp-1, data=dat, constr=cm1, rhs=0, nec=0)
#' 
#' # order restriction (increasing means)
#' cm2 <- rbind(c(-1,1,0,0),
#'              c(0,-1,1,0),
#'              c(0,0,-1,1))
#' fm2 <- orlm(y ~ grp-1, data=dat, constr=cm2,
#'             rhs=rep(0,nrow(cm2)), nec=0)
#' 
#' # order restriction (increasing at least by delta=1)
#' fm3 <- orlm(y ~ grp-1, data=dat, constr=cm2,
#'             rhs=rep(1,nrow(cm2)), nec=0)
#'
#' # larger than average of the neighboring first 2 parameters
#' cm4 <- rbind(c(-0.5,-0.5,1,0),
#'              c(0,-0.5,-0.5,1))
#' fm4 <- orlm(y ~ grp-1, data=dat, constr=cm4,
#'             rhs=rep(0,nrow(cm4)), nec=0)
#' 
#' # equality constraints (all parameters equal)
#' fm5 <- orlm(y ~ grp-1, data=dat, constr=cm2,
#'             rhs=rep(0,nrow(cm2)), nec=nrow(cm2))
#'             
#' # alternatively
#' fm5 <- orlm(y ~ grp-1, data=dat, constr=cm2,
#'             rhs=rep(0,nrow(cm2)), nec=c(TRUE,TRUE,TRUE))
#'             
#' # constraining the 1st and the 4th parameter
#' # to their true values, and the 2nd and 3rd between them
#' cm6 <- rbind(c( 1,0,0,0),
#'              c(-1,1,0,0),
#'              c(0,-1,0,1),
#'              c(-1,0,1,0),
#'              c(0,0,-1,1),
#'              c(0,0, 0,1))
#' fm6 <- orlm(y ~ grp-1, data=dat, constr=cm6,
#'             rhs=c(1,rep(0,4),5), nec=c(TRUE,rep(FALSE,4),TRUE))
#'             
#'             
#' ###############################################################
#' ## Example from Kuiper, R.M. and Hoijtink, H. (Unpublished). ##
#' ## A Fortran 90 program for the generalization of the        ##
#' ## order restricted information criterion.                   ##
#' ###############################################################
#' 
#' # constraint definition
#' cmat <- cbind(diag(3), 0) + cbind(0, -diag(3))
#' constr <- kronecker(diag(3), cmat)
#' 
#' # no effect model
#' (fm0 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
#'              constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr)))
#' 
#' # order constrained model (increasing serum levels with increasing doses)
#' fm1 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
#'             constr=constr, rhs=rep(0, nrow(constr)), nec=0)
#' summary(fm1)
#' 
#' # unconstrained model
#' (fmunc <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
#'                constr=matrix(0, nrow=1, ncol=12), rhs=0, nec=0))


orlm <-
function(formula, data, constr, rhs, nec, control=orlmcontrol()){
UseMethod("orlm")
}

#' @rdname orlm
orlm.formula <-
function(formula, data, constr, rhs, nec, control=orlmcontrol()){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- cbind(model.response(mf, "numeric"))
  x <- model.matrix(mt, mf, contrasts)

  if (is.numeric(constr)) constr <- rbind(constr)
  if (!is.matrix(constr)) stop("constr needs to be a matrix.")
  if (ncol(y)*ncol(x) != ncol(constr)) stop(paste("constr has not correct dimensions.\nNumber of columns (",ncol(constr),") should equal the number of responses times the number of parameters: ",ncol(y), "*", ncol(x), "=",ncol(y)*ncol(x), sep=""))
  if (length(rhs) != nrow(constr)) stop(paste("rhs has a different number of elements than there are numbers of rows in constr (",length(rhs), " != ", nrow(constr), ")", sep=""))
  if (is.numeric(nec) & length(nec) != 1) stop("nec needs to be single a numeric value or a logical vector with the same length as the number of constraints.")
  if (is.logical(nec) & length(nec) != length(rhs)) stop("nec needs to be single a numeric value or a logical vector with the same length as the number of constraints.")
  if (is.logical(nec)){
    ord <- order(nec, decreasing=TRUE)
    constr <- constr[ord,,drop=FALSE]
    rhs <- rhs[ord]
    nec <- sum(nec)
  }
  if (nec < 0) stop("nec needs to be positive")
  if (nec > length(rhs)) stop(paste("nec is larger than the number of constraints. (",nec," > ",length(rhs),")", sep=""))    

  ########################
  ## unconstrained linear model
  unc <- lm(y ~ x-1)

  tBeta <- as.vector(coefficients(unc))
  Sigma <- (t(residuals(unc)) %*% residuals(unc))/nrow(x)
  detU <- 1

  ############################
  ## lin model with order restrictions

  orsolve <- function(tBeta, x, y, Constr, RHS, NEC){
    Sigma <- (t(y - x %*% matrix(tBeta, ncol=ncol(y))) %*% (y - x %*% matrix(tBeta, ncol=ncol(y))))/nrow(x)
    yVx <- kronecker(solve(Sigma), t(x)) %*% as.vector(y)
    dvec <- 2*yVx
    Dmat <- 2*kronecker(solve(Sigma), t(x) %*% x)
    Amat <- t(Constr)
    bvec <- RHS
    solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=NEC)    
  }

  orBeta <- tBeta
  val <- 0
  for (i in 1:control$maxiter){
    sqp <- orsolve(orBeta, x, y, constr, rhs, nec)
    orBeta <- sqp$solution
    if (abs(sqp$value - val) <= control$absval) break else val <- sqp$value
  }
  if (i == control$maxiter & abs(sqp$value - val) > control$absval) warning("Maximum number of iterations reached without convergence!")
  orSigma <- (t(y - x %*% matrix(orBeta, ncol=ncol(y))) %*% (y - x %*% matrix(orBeta, ncol=ncol(y))))/nrow(x)

  N <- length(as.vector(y))
  loglik <- (-N/2.0)*log(2*pi) + (-1/2.0)*(nrow(x)*log(det(orSigma)) + ncol(y)*log(detU)) - (1/2.0)*N

  if (ncol(y) > 1){
    orBeta <- matrix(orBeta, ncol=ncol(y), dimnames=list(colnames(x), colnames(y)))
    tBeta <- matrix(tBeta, ncol=ncol(y), dimnames=list(colnames(x), colnames(y)))
  }
  names(orBeta) <- colnames(x)
  out <- list(call=cl, X=x, y=y, unccoefficients=tBeta, coefficients=orBeta, fitted=x %*% orBeta, residuals=y - x %*% orBeta, sigma=Sigma, orSigma=orSigma, logLik=loglik, constr=constr, rhs=rhs, nec=nec, Niter=i, iact=sqp$iact)
  class(out) <- "orlm"
  return(out)
}

