\name{goric}
\alias{goric}
\alias{goric.orlm}
\alias{goric.orgls}
\alias{goric.orglm}
\alias{goric.list}
\alias{goric_penalty}
\alias{orglm_penalty}

\title{Calculate GORIC}
\description{
It calculates the order-restricted log likelihood, the penalty of the GORIC,
the GORIC values, differences to the minimum GORIC value, and the GORIC
weights for a set of hypotheses (orlm objects), where the penalty is based on \eqn{iter} iterations.
}
\usage{
goric(object, ..., iter = 100000, type="GORIC", dispersion=1, mc.cores=1)
\method{goric}{orlm}(object, ..., iter = 100000, type="GORIC", mc.cores=1)
\method{goric}{orgls}(object, ..., iter = 100000, type="GORIC", mc.cores=1)
\method{goric}{list}(object, ..., iter = 100000, type="GORIC", dispersion=1, mc.cores=1)
\method{goric}{orglm}(object, ..., iter = 100000, type="GORIC", dispersion=1, mc.cores=1)
goric_penalty(object, iter = 1e+05, type="GORIC", mc.cores=1)
orglm_penalty(object, iter = 100000, type="GORIC", mc.cores=1)
}
\arguments{
  \item{object}{an object of class orlm, orgls, orglm, or a list of these objects}
  \item{\dots}{further objects of class orlm, orgls, OR orglm}
  \item{iter}{number of iterations to calculate GORIC penalty terms}
  \item{type}{if \code{"GORIC"} (default), the penalty term for the
    generalized order restriction is computed; with \code{"GORICCa"} or
    \code{"GORICCb"} small sample corrections for the penalty term are
    applied (only available for orlm models)}
  \item{dispersion}{dispersion parameter to scale GORIC analogue to QAIC
  in generalized linear models}
  \item{mc.cores}{number of cores using a socket cluster implemented in package
  \code{parallel}}
}
\details{
The hypothesis with the lowest GORIC value is the preferred one.

The GORIC weights reflect the support of each hypothesis in the set.
To compare two hypotheses (and not one to the whole set), one should examine 
the ratio of the two corresponding GORIC weights.


To safequard for weak hypotheses (i.e., hypotheses not supported by the data), 
one should include a model with no constraints (the so-called unconstrained model).
More information can be obtained from Kuiper, Hoijtink, and Silvapulle (2011) and 
Kuiper, Hoijtink, and Silvapulle (unpublished).
}
\value{a data.frame with the information criteria or a single penalty term}
\references{
Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
\bold{98}, 495--501.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the Order-Restricted Information Criterion for Multivariate Normal Linear Models. \emph{Journal of Statistical Planning and Inference}, \bold{142}, 2454-2463. doi:10.1016/j.jspi.2012.03.007.

Kuiper R.M. and Hoijtink H. (submitted). A Fortran 90 Program for the Generalization of the Order-Restricted Information Criterion. Journal of Statictical Software.
}
\author{
Daniel Gerhard and Rebecca M. Kuiper
}


\seealso{\code{\link{orlm}}, \code{\link{orgls}}}
\examples{
## Example from Kuiper, R.M. and Hoijtink, H. (Unpublished).
# A Fortran 90 program for the generalization of the
# order restricted information criterion.

# constraint definition
cmat <- cbind(diag(3), 0) + cbind(0, -diag(3))
constr <- kronecker(diag(3), cmat)
constr

# no effect model
fm0 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
            constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr))
fm0

# order constrained model (increasing serum levels with increasing doses)
fm1 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
            constr=constr, rhs=rep(0, nrow(constr)), nec=0)
summary(fm1)

# unconstrained model
fmunc <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
              constr=matrix(0, nrow=1, ncol=12), rhs=0, nec=0)
fmunc

# calculate GORIC
# (only small number of iterations to decrease computation time, default: iter=100000)
goric(fm0, fm1, fmunc, iter=1000)
}

\keyword{models}
