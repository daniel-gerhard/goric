#' Simulation from order restricted linear models
#' 
#' Simulation function for orlm and orgls objects
#' 
#' Given the estimated coefficients of a orlm or orgls model, a set new parameters are generated. n.sims new sets of observations are generated based on the unrestricted model; these new datasets are used to estimate a new set of model coefficients incorporating the given order restrictions.
#' 
#' @param object an object of class "orlm" or "orgls".
#' @param n.sims number of simulation replications.
#' 
#' @return a list with sets of simulated parameters.
#' 
#' @seealso \code{\link{orlm}}, \code{\link{orgls}}
#' 
#' @keywords methods
#' 
#' @examples 
#' ########################
#' ## Artificial example ##
#' ########################
#' n <- 10
#' m <- c(1,1,2)
#' dat <- data.frame(grp=as.factor(rep(1:length(m), each=n)),
#'                   y=rnorm(n*length(m), rep(m, each=n), 1))
#' cm <- rbind(c(-1,1,0),
#'             c(0,-1,1))
#' fm <- orlm(y ~ grp-1, data=dat, constr=cm, rhs=rep(0,nrow(cm)), nec=0)
#' b <- sim(fm, n.sims=1000)$coef
#' pairs(t(b), cex=0.3)

sim <- function(object, n.sims){
  UseMethod("sim")
}

#' @rdname sim
sim.orlm <- function(object, n.sims){
  object <- object
  unr <- lm(object$y ~ object$X-1)
  ### code taken from package arm
  summ <- summary(unr)
  coef <- summ$coef[,1:2,drop=FALSE]
  dimnames(coef)[[2]] <- c("coef.est","coef.sd")
  sigma.hat <- summ$sigma
  beta.hat <- coef[,1,drop = FALSE]
  V.beta <- summ$cov.unscaled
  n <- summ$df[1] + summ$df[2]
  k <- summ$df[1]
  sigma <- sigma.hat*sqrt((n-k)/rchisq(n.sims,n-k))
  # new simulated dataset
  smat <- matrix(rep(object$X %*% coefficients(unr), n.sims), ncol=n.sims) + sapply(sigma, function(sigma) rnorm(nrow(object$X), 0, sigma))
  # new set of constraint parameters
  beta <- rbind(apply(smat, 2, function(sy) coefficients(orlm(sy ~ object$X-1, constr=object$constr, rhs=object$rhs, nec=object$nec))))
  rownames(beta) <- colnames(object$X)
  return(list(coef=beta, sigma=sigma))
}

#' @rdname sim
sim.orgls <- function(object, n.sims){
  object <- object
  unr <- lm(object$y ~ object$X-1)
  ### code taken from package arm
  summ <- summary(unr)
  coef <- summ$coef[,1:2,drop=FALSE]
  dimnames(coef)[[2]] <- c("coef.est","coef.sd")
  sigma.hat <- summ$sigma
  beta.hat <- coef[,1,drop = FALSE]
  V.beta <- summ$cov.unscaled
  n <- summ$df[1] + summ$df[2]
  k <- summ$df[1]
  sigma <- sigma.hat*sqrt((n-k)/rchisq(n.sims,n-k))
  # new simulated dataset
  smat <- matrix(rep(object$X %*% coefficients(unr), n.sims), ncol=n.sims) + sapply(sigma, function(sigma) rnorm(nrow(object$X), 0, sigma))
  # new set of constraint parameters
  beta <- rbind(apply(smat, 2, function(sy) coefficients(orlm(sy ~ object$X-1, constr=object$constr, rhs=object$rhs, nec=object$nec))))
  rownames(beta) <- colnames(object$X)
  return(list(coef=beta, sigma=sigma))
}
