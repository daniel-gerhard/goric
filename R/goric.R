#' Calculate GORIC
#' 
#' The goric function calculates the order-restricted log likelihood, the penalty of the generalised order restricted information criterion (GORIC), the GORIC values, differences to the minimum GORIC value, and the GORIC weights for a set of hypotheses, where the penalty is based on \eqn{iter} iterations.
#' The hypothesis with the lowest GORIC value is the preferred one. 
#' The GORIC weights reflect the support of each hypothesis in the set. To compare two hypotheses (and not one to the whole set), one should examine the ratio of the two corresponding GORIC weights.
#' To safequard for weak hypotheses (i.e., hypotheses not supported by the data), one should include a model with no constraints (the so-called unconstrained model).
#' 
#' @param object an object of class orlm, orgls, orglm, or a list of these objects
#' @param ... further objects of class orlm, orgls, or orglm
#' @param iter number of iterations to calculate GORIC penalty terms
#' @param type if \code{"GORIC"} (default), the penalty term for the generalized order restriction information criterion is computed; with \code{"GORICCa"} or \code{"GORICCb"} small sample corrections for the penalty term are applied
#' @param dispersion dispersion parameter to scale GORIC analogously to QAIC in generalized linear models
#' @param mc.cores number of cores using a socket cluster implemented in package \code{parallel}
#' 
#' @return a data.frame with the information criteria or a single penalty term
#' 
#' @references 
#' \itemize{
#' \item Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, \bold{98}, 495--501.
#' \item Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the Order-Restricted Information Criterion for Multivariate Normal Linear Models. \emph{Journal of Statistical Planning and Inference}, \bold{142}, 2454-2463. doi:10.1016/j.jspi.2012.03.007.
#' \item Kuiper R.M. and Hoijtink H. (submitted). A Fortran 90 Program for the Generalization of the Order-Restricted Information Criterion. Journal of Statictical Software.
#' }
#' 
#' @seealso \code{\link{orlm}}, \code{\link{orgls}}
#' 
#' @keywords models
#' 
#' @examples 
#' ## Example from Kuiper, R.M. and Hoijtink, H. (Unpublished).
#' # A Fortran 90 program for the generalization of the 
#' # order restricted information criterion.
#' # constraint definition
#' cmat <- cbind(diag(3), 0) + cbind(0, -diag(3))
#' constr <- kronecker(diag(3), cmat)
#' constr
#' 
#' # no effect model
#' (fm0 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
#'             constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr)))
#' 
#' # order constrained model (increasing serum levels with increasing doses)
#' fm1 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
#'             constr=constr, rhs=rep(0, nrow(constr)), nec=0)
#' summary(fm1)
#' 
#' # unconstrained model
#' (fmunc <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
#'               constr=matrix(0, nrow=1, ncol=12), rhs=0, nec=0))
#' 
#' # calculate GORIC
#' # (only small number of iterations to decrease computation time, default: iter=100000)
#' goric(fm0, fm1, fmunc, iter=1000)

goric <-
function(object, ..., iter=100000, type="GORIC", dispersion=1, mc.cores=1){
  UseMethod("goric")
}

#' @rdname goric
goric.orlm <-
function(object, ..., iter=100000, type="GORIC", mc.cores=1){
  if (!inherits(object, "orlm") & !inherits(object, "list")) stop("object needs to be of class orlm or a list of orlm objects")
  if (iter < 1) stop("No of iterations < 1")
  if (inherits(object, "orlm")) objlist <- list(object, ...) else objlist <- object
  isorlm <- sapply(objlist, function(x) inherits(x, "orlm"))
  orlmlist <- objlist[isorlm]  
  Call <- match.call()
  Call$iter <- NULL
  Call$type <- NULL
  Call$mc.cores <- NULL
  if (inherits(object, "orlm")) names(orlmlist) <- as.character(Call[-1L])[isorlm]
  loglik <- sapply(orlmlist, function(x) x$logLik)
  penalty <- sapply(orlmlist, function(x) goric_penalty(x, iter=iter, type=type, mc.cores=mc.cores))
  goric <- -2*(loglik - penalty)
  delta <- goric - min(goric)
  goric_weights <- exp(-delta/2) / sum(exp(-delta/2))
  data.frame(loglik, penalty, goric=goric, goric_weights=round(goric_weights,3))
}

#' @rdname goric
goric.orgls <-
function(object, ..., iter=100000, type="GORIC", mc.cores=1){
  if (!inherits(object, "orgls") & !inherits(object, "list")) stop("object needs to be of class orgls or a list of orgls objects")
  if (type != "GORIC") stop("Only type='GORIC' is implemented for orgls objects!")
  if (iter < 1) stop("No of iterations < 1")
  if (inherits(object, "orgls")) objlist <- list(object, ...) else objlist <- object
  isorgls <- sapply(objlist, function(x) inherits(x, "orgls"))
  orlmlist <- objlist[isorgls]  
  Call <- match.call()
  Call$iter <- NULL
  Call$type <- NULL
  Call$mc.cores <- NULL
  if (inherits(object, "orgls")) names(orlmlist) <- as.character(Call[-1L])[isorgls]
  loglik <- sapply(orlmlist, function(x) x$logLik)
  ep <- sapply(orlmlist, function(x) x$extrap)
  penalty <- sapply(orlmlist, function(x) goric_penalty(x, iter=iter, type=type, mc.cores=mc.cores))
  goric <- -2*(loglik - penalty - ep)
  delta <- goric - min(goric)
  goric_weights <- exp(-delta/2) / sum(exp(-delta/2))
  data.frame(loglik, penalty, vcdf=ep, goric=goric, goric_weights=round(goric_weights,3))
}

#' @rdname goric
goric.list <- function(object, ..., iter=100000, type="GORIC", dispersion=1, mc.cores=1){
  if (all(sapply(object, class) == "orlm")) out <- goric.orlm(object, iter=iter, type=type, mc.cores=mc.cores)
  if (all(sapply(object, class) == "orgls")) out <- goric.orgls(object, iter=iter, type=type, mc.cores=mc.cores)
  if (all(sapply(object, class) == "orglm")) out <- goric.orglm(object, iter=iter, type=type, dispersion=dispersion, mc.cores=mc.cores)
  return(out)
}

#' @rdname goric
goric.orglm <- function(object, ..., iter=100000, type="GORIC", dispersion=1, mc.cores=1){
  if (inherits(object, "orglm")) objlist <- list(object, ...) else objlist <- object
  isorglm <- sapply(objlist, function(x) inherits(x, "orglm"))
  orglmlist <- objlist[isorglm]
  Call <- match.call()
  Call$iter <- NULL
  Call$type <- NULL
  Call$dispersion <- NULL
  Call$mc.cores <- NULL
  if (inherits(object, "orglm")) names(orglmlist) <- as.character(Call[-1L])[isorglm]
  loglik <- sapply(orglmlist, function(x){
    fam <- x$family
    p <- x$rank
    if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) p <- p + 1
    return(p - x$oaic/2)
  })  
  penalty <- sapply(orglmlist, function(x){
    penal <- orglm_penalty(object=x, iter=iter, type=type, mc.cores=mc.cores)
    fam <- x$family
    if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) penal <- penal + 1
    return(penal)
  })
  
  goric <- -2*loglik/dispersion + 2*penalty
  delta <- goric - min(goric)
  goric_weights <- exp(-delta/2)/sum(exp(-delta/2))
  data.frame(loglik, penalty, goric = goric, goric_weights = round(goric_weights, 3))
}

