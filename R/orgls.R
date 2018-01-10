#' Fitting generalized least squares regression models with order restrictions
#' 
#' \code{orgls} is used to fit generalised least square models analogously to the function \code{gls} in package \code{nlme} but with order restrictions on the parameters.
#' 
#' The contraints in the hypothesis of interest are defined by \eqn{constr}, \eqn{rhs}, and \eqn{nec}. The first \eqn{nec} constraints are the equality contraints: \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}; and the remaing ones are the inequality contraints: \eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}.
#' Two requirements should be met:
#' \enumerate{
#'     \item The first \eqn{nec} constraints must be the equality contraints (i.e., \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}) and the remaining ones the inequality contraints (i.e., \eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}).
#'     \item When \eqn{rhs} is not zero, \eqn{Constr} should be of full rank (after discarding redundant restrictions).
#' }
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which orgls is called.
#' @param constr matrix with constraints; with rows as constraint definition, columns should be in line with the parameters of the model
#' @param rhs vector of right hand side elements; \eqn{Constr \; \theta \geq rhs}; number should equal the number of rows of the constr matrix
#' @param nec number of equality constraints; a numeric value treating the first nec constr rows as equality constraints, or a logical vector with \code{TRUE} for equality- and \code{FALSE} for inequality constraints.
#' @param weights a \code{\link{varClasses}} object; more details are provided on the help pages in R package \code{nlme}
#' @param correlation a \code{\link{corClasses}} object; more details are provided on the help pages in R package \code{nlme}
#' @param control a list of control arguments; see \code{\link{orlmcontrol}} for details.
#' 
#' @return an object of class orgls
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
#' # generating example data
#' library(mvtnorm)
#' # group means
#' m <- c(0,5,5,7)
#' # compound symmetry structure of residuals
#' # (10 individuals per group, rho=0.7) 
#' cormat <- kronecker(diag(length(m)*10), matrix(0.7, nrow=length(m), ncol=length(m)))
#' diag(cormat) <- 1
#' # different variances per group
#' sds <- rep(c(1,2,0.5,1), times=10*length(m))
#' sigma <- crossprod(diag(sds), crossprod(cormat, diag(sds)))
#' response <- as.vector(rmvnorm(1, rep(m, times=10*length(m)), sigma=sigma))
#' dat <- data.frame(response,
#'                   grp=rep(LETTERS[1:length(m)], times=10*length(m)), 
#'                   ID=as.factor(rep(1:(10*length(m)), each=length(m))))
#'                   
#' ## set of gls models:
#' # unconstrained model
#' m1 <- orgls(response ~ grp-1, data = dat,
#'             constr=rbind(c(0,0,0,0)), rhs=0, nec=0,
#'             weights=varIdent(form=~1|grp),
#'             correlation=corCompSymm(form=~1|ID))
#'
#' # simple order
#' m2 <- orgls(response ~ grp-1, data = dat,
#'             constr=rbind(c(-1,1,0,0),c(0,-1,1,0),c(0,0,-1,1)), rhs=c(0,0,0), nec=0,
#'             weights=varIdent(form=~1|grp),
#'             correlation=corCompSymm(form=~1|ID))
#' 
#' # equality constraints
#' m3 <- orgls(response ~ grp-1, data = dat,
#'             constr=rbind(c(-1,1,0,0),c(0,-1,1,0),c(0,0,-1,1)), rhs=c(0,0,0), nec=3,
#'             weights=varIdent(form=~1|grp),
#'             correlation=corCompSymm(form=~1|ID))


orgls <-
function(formula, data, constr, rhs, nec, weights=NULL, correlation=NULL, control=orlmcontrol()){
UseMethod("orgls")
}

#' @rdname orgls
orgls.formula <-
function(formula, data, constr, rhs, nec, weights=NULL, correlation=NULL, control=orlmcontrol()){
  cl <- match.call()
  if (!is.null(correlation)){
    groups <- getGroupsFormula(correlation)
  } else {
    groups <- NULL
  }
  glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
  model <- terms(formula, data = data)
  mfArgs <- list(formula = asOneFormula(formula(glsSt), formula, groups), data = data, na.action = na.fail)
  dataMod <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMod)
  if (!is.null(groups)) {
    groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), sep = "|")))
    grps <- getGroups(dataMod, groups, level = length(getGroupsFormula(groups, asList = TRUE)))
    ord <- order(grps)
    grps <- grps[ord]
    dataMod <- dataMod[ord, , drop = FALSE]
    revOrder <- match(origOrder, row.names(dataMod))
  } else {
    grps <- NULL
  }
  X <- model.frame(model, dataMod)
  contr <- lapply(X, function(el) if (inherits(el, "factor")) contrasts(el))
  contr <- contr[!unlist(lapply(contr, is.null))]
  x <- model.matrix(model, X)
  y <- eval(model[[2]], dataMod)

  if (is.numeric(constr)) constr <- rbind(constr)
  if (!is.matrix(constr)) stop("constr needs to be a matrix.")
  if (ncol(x) != ncol(constr)) stop(paste("constr has not correct dimensions.\nNumber of columns (",ncol(constr),") should equal the number of parameters: ", ncol(x), sep=""))
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
  unc <- gls(formula, data=dataMod, weights=weights, correlation=correlation, method="ML")

  ## extracting the variance-covariance structure
  if (is.null(unc$modelStruct$varStruct)){
    V <- diag(nrow(x))
  } else {
    V <- diag(attr(unc$modelStruct$varStruct, "weights"))
  }
  if (is.null(unc$modelStruct$corStruct)){
    crr <- diag(nrow(x))
  } else {
    cr <- corMatrix(unc$modelStruct$corStruct)
    crr <- if (is.matrix(cr)) cr else as.matrix(bdiag(cr))
  }
  W <- V %*% crr %*% V

  tBeta <- lm.gls(formula, data = dataMod, W=W)$coefficients
  
  # taken from lm.gls in package MASS
  # transforming X and y into a classical linear model framework
  eW <- eigen(W, TRUE)
  d <- eW$values
  if (any(d <= 0)) stop("'W' is not positive definite")
  eWv <- eW$vector
  A <- diag(sqrt(d)) %*% t(eWv)
  Ainv <- eWv %*% diag(1/sqrt(d))
  X <- A %*% x
  Y <- as.vector(A %*% y)  
  res <- Y - X %*% tBeta
  Sigma <- as.vector(t(res) %*% (res))/(nrow(x))  
  ############################
  ## lin model with order restrictions
  orsolve <- function(tBeta, X, Y, Constr, RHS, NEC){
    yVx <- t(X) %*% Y
    dvec <- 2*yVx
    Dmat <- 2*(t(X) %*% X)
    Amat <- t(Constr)
    bvec <- RHS
    solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=NEC)    
  }
  orBeta <- tBeta
  val <- 0
  for (i in 1:control$maxiter){
    sqp <- orsolve(orBeta, X, Y, constr, rhs, nec)
    orBeta <- sqp$solution
    if (abs(sqp$value - val) <= control$absval) break else val <- sqp$value
  }
  if (i == control$maxiter & abs(sqp$value - val) > control$absval) warning("Maximum number of iterations reached without convergence!")

  ores <- (Y - X %*% orBeta)
  orSigma <- as.vector(t(ores) %*% (ores))/(nrow(x))
  Aores <- Ainv %*% (Y - X %*% orBeta)
  AorSigma <- as.vector(t(Aores) %*% diag(diag(W)) %*% (Aores))/(nrow(x))  
  p <- unc$dims$p
  N <- unc$dims$N
  Np <- N - p
  loglik <- (-N/2)*log(2*pi) + (-1/2)*(nrow(x)*log(AorSigma) + determinant(W, logarithm=TRUE)$modulus) - (1/2)*N + sum(log(diag(W)))

  names(orBeta) <- colnames(x)
  out <- list(call=cl, X=X, XW=x, y=Y, unccoefficients=tBeta, coefficients=orBeta, fitted=Ainv %*% (X %*% orBeta), residuals=Ainv %*% (Y - X %*% orBeta), sigma=Sigma, orSigma=orSigma, logLik=loglik, constr=constr, rhs=rhs, nec=nec, Niter=i, iact=sqp$iact, extrap=length(coef(unc[["modelStruct"]])), modelStruct=unc$modelStruct, W=W)
  class(out) <- "orgls"
  return(out)
}
