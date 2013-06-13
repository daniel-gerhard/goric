\name{orlmSet}
\alias{orlmSet}

\title{Fitting multivariate regression models with order restrictions}
\description{It determines the order-restricted maximum likelihood estimates 
and the corresponding log likelihood for the hypothesis of interest.
Additionally it gives the (unconstrained) maximum likelihood estimates and
the active contraints.}

\usage{
orlmSet(formula, data, set, direction="increase", n=NULL, base=1,
        control = orlmcontrol())
}
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
}
  \item{set}{either a character string (see \code{\link{constrSet}}), or
  a list with slots for constr, rhs, and nec similarly defined as in
  \code{\link{orlm}}}
\item{direction}{direction of the order constraints}
\item{n}{a (possibly named) vector of sample sizes for each group}
\item{base}{column of the constraint matrix representing a control
    group}
  \item{control}{a list of control arguments; see
    \code{\link{orlmcontrol}} for details.}
}
\details{
  This function is just a wrapper for repeated calls of
  \code{\link{orlm}} with different constraint definitions.
  Predefined lists with constraint-sets can be constructed with function \code{\link{constrSet}}.
}
\value{an object of class list}
\references{
Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
\bold{98}, 495--501.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the Order-Restricted Information Criterion for Multivariate Normal Linear Models. \emph{Journal of Statistical Planning and Inference}, \bold{142}, 2454-2463. doi:10.1016/j.jspi.2012.03.007.

Kuiper R.M. and Hoijtink H. (submitted). A Fortran 90 Program for the Generalization of the Order-Restricted Information Criterion. Journal of Statictical Software.
}
\author{Daniel Gerhard and Rebecca M. Kuiper}


\seealso{\code{\link{orlm}}, \code{\link{constrSet}}, \code{\link{goric}}}

\examples{
########################
## Artificial example ##
########################

n <- 10
m <- c(1,2,4,5,2,1)
nm <- length(m)
dat <- data.frame(grp=as.factor(rep(1:nm, each=n)),
                  y=rnorm(n*nm, rep(m, each=n), 1))

(cs <- constrSet(table(dat$grp), set="sequence"))
oss <- orlmSet(y ~ grp-1, data=dat, set=cs)
oss
# the same as:
oss <- orlmSet(y ~ grp-1, data=dat, set="sequence")
}

\keyword{models}
