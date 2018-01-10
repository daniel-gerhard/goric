#' Generate Constraint Sets
#' 
#' Generate sets of constraint matrices (constr), right hand side elements, and numbers of equality constraints (nec) with a predefined structure
#' 
#' @param n a (possibly named) vector of sample sizes for each group.
#' @param set character string defining the type of constraints; one of "sequence", "seqcontrol", "lplateau", "uplateau", or "downturn"
#' @param direction direction of the inequality constraints, either "increase" or "decrease"
#' @param base column of the constraint matrix representing a control group
#' 
#' @return a list with slots constr, rhs, and nec for each constraint definition
#' 
#' @seealso \code{\link{orlm}}, \code{\link{constrMat}}
#' 
#' @keywords misc
#' 
#' @examples 
#' n <- c(10,20,30,40)
#' constrSet(n, set="sequence")
#' constrSet(n, set="seqcontrol")
#' constrSet(n, set="lplateau")
#' constrSet(n, set="uplateau")
#' constrSet(n, set="downturn")
#' constrSet(n, set="williams")

constrSet <- function(n, set=c("sequence","seqcontrol","lplateau","uplateau","downturn","williams"), direction=c("increase","decrease"), base=1){
  d <- length(n)
  if (d < 2) stop("less than two groups")
  if (!is.numeric(n)) stop(sQuote("n"), " is not numeric")
  if (base < 1 || base > d) stop("base is not between 1 and ", d)
  if (!is.null(names(n))) varnames <- names(n) else varnames <- 1:d
  direction <- match.arg(direction)
  switch(direction, increase = {
    direct <- 1
    vorz <- "<"
  }, decrease = {
    direct <- -1
    vorz <- ">"
  })
  set <- match.arg(set)
  switch(set, sequence = {
    sl <- lapply(1:(d-1), function(i){
      list(constr=direct*constrMat(n, type="monotone")[1:i,,drop=FALSE],
           rhs=rep(0, i),
           nec=rep(FALSE,i) )
    })
    CM <- matrix(0,nrow=1, ncol=d)
    colnames(CM) <- varnames
    sl[[d]] <- list(constr=CM, rhs=0, nec=FALSE)
    nams <- sapply(2:d, function(i) paste(varnames[1:i], collapse=vorz))
    names(sl) <- c(nams, "unconstrained")
  }, seqcontrol = {
    sl <- lapply(1:(d-1), function(i){
      list(constr=direct*constrMat(n, type="control", base=base)[1:i,,drop=FALSE],
           rhs=rep(0, i),
           nec=rep(FALSE,i) )
    })
    sl[[d]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply(1:(d-1), function(i) paste(paste(varnames[base], varnames[(1:d)[-base]][1:i], sep=vorz), collapse="|"))
    names(sl) <- c(nams, "unconstrained")
  }, lplateau = {
    sl <- lapply(1:d, function(i){
      list(constr=direct*constrMat(n, type="monotone"),
           rhs=rep(0, d-1),
           nec=c(rep(TRUE,d-i), rep(FALSE,i-1)) )
    })
    sl[[d+1]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply(1:d, function(i){
      vz <- rep(vorz, length=d-1)
      vz[sl[[i]]$nec] <- "="
      vz <- c(vz, "")
      paste(paste(varnames[1:d], vz, sep=""), collapse="")
    })
    names(sl) <- c(nams, "unconstrained")
  }, uplateau = {
    sl <- lapply(1:d, function(i){
      list(constr=direct*constrMat(n, type="monotone"),
           rhs=rep(0, d-1),
           nec=c(rep(FALSE,d-i), rep(TRUE,i-1)) )
    })
    sl[[d+1]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply(1:d, function(i){
      vz <- rep(vorz, length=d-1)
      vz[sl[[i]]$nec] <- "="
      vz <- c(vz, "")
      paste(paste(varnames[1:d], vz, sep=""), collapse="")
    })
    names(sl) <- c(nams, "unconstrained")
  }, downturn = {
    sl <- list()
    sl[2:(d-1)] <- lapply((d-1):2, function(i){
      cm <- constrMat(n, type="monotone")
      cm[i:(d-1),] <- -1*cm[i:(d-1),]
      list(constr=cm,
           rhs=rep(0, d-1),
           nec=rep(FALSE,d-1) )
    })
    sl[[1]] <- list(constr=direct*constrMat(n, type="monotone"), rhs=rep(0, d-1),nec=rep(FALSE,d-1))
    sl[[d]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply((d-1):1, function(i){
      if (vorz == ">") vz <- rep("<", length=d-1) else vz <- rep(">", length=d-1)
      vz[1:i] <- vorz
      vz <- c(vz, "")
      paste(paste(varnames[1:d], vz, sep=""), collapse="")
    })
    names(sl) <- c(nams, "unconstrained")
  }, williams = {
    sl <- lapply(1:(d-1), function(i){
      list(constr=direct*constrMat(n, type="caverage", base=1)[1:i,,drop=FALSE],
           rhs=rep(0, i),
           nec=rep(FALSE,i) )
    })
    CM <- matrix(0,nrow=1, ncol=d)
    colnames(CM) <- varnames
    sl[[d]] <- list(constr=CM, rhs=0, nec=FALSE)
    nams <- sapply(2:d, function(i) paste(varnames[1], " ", vorz, " ave(", paste(varnames[2:i], collapse=","), ")", sep=""))
    nams[1] <- paste(varnames[1], vorz, varnames[2], sep=" ")
    names(sl) <- c(nams, "unconstrained")
  })
  return(sl)
}
