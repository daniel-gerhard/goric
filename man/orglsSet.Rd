\name{orglsSet}
\alias{orglsSet}

\title{Fitting generalized least squares regression models with order restrictions}
\description{It determines the order-restricted maximum likelihood estimates 
and the corresponding log likelihood for the hypothesis of interest.
Additionally it gives the (unconstrained) maximum likelihood estimates and
the active contraints.}

\usage{
orglsSet(formula, data, weights=NULL,
         correlation=NULL, set, direction="increase",
         n=NULL, base=1, control = orlmcontrol())
}
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
  }
    \item{weights}{a \code{\link{varClasses}} object; more details are
    provided on the help pages in R package \code{nlme}}
  \item{correlation}{a \code{\link{corClasses}} object; more details are
    provided on the help pages in R package \code{nlme}}
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
  \code{\link{orgls}} with different constraint definitions.
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


\keyword{models}
