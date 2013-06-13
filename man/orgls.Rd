\name{orgls}
\alias{orgls}
\alias{orgls.formula}

\title{Fitting generalized least squares regression models with order restrictions}
\description{It determines the order-restricted maximum likelihood estimates 
and the corresponding log likelihood for the hypothesis of interest.
Additionally it gives the (unconstrained) maximum likelihood estimates and
the active contraints.}

\usage{
orgls(formula, data, constr, rhs, nec, weights=NULL, correlation=NULL,
      control = orlmcontrol())
\method{orgls}{formula}(formula, data, constr, rhs, nec, weights=NULL,
correlation=NULL,
      control = orlmcontrol())
}
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which orgls is called.
}
  \item{constr}{matrix with constraints; with rows as constraint
    definition, columns should be in line with the parameters of the model}
  \item{rhs}{vector of right hand side elements; \eqn{Constr \; \theta \geq
      rhs}; number should equal the number of rows of the constr matrix}
  \item{nec}{number of equality constraints; a numeric value treating
    the first nec constr rows as equality constraints, or a logical
    vector with \code{TRUE} for equality- and \code{FALSE} for
    inequality constraints.}
  \item{weights}{a \code{\link{varClasses}} object; more details are
    provided on the help pages in R package \code{nlme}}
  \item{correlation}{a \code{\link{corClasses}} object; more details are
    provided on the help pages in R package \code{nlme}}
  \item{control}{a list of control arguments; see \code{\link{orlmcontrol}} for details.
  }
}
\details{
The contraints in the hypothesis of interest are defined by 
\eqn{Constr}, \eqn{rhs}, and \eqn{nec}.
The first \eqn{nec} constraints are the equality contraints:
\eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]};
and the remaing ones are the inequality contraints:
\eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}.

Two requirements should be met:
\enumerate{
\item The first \eqn{nec} constraints must be the equality contraints 
(i.e., \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}) 
and the remaining ones the inequality contraints
(i.e., \eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}).
\item When \eqn{rhs} is not zero, \eqn{Constr} should be of full rank 
(after discarding redundant restrictions). 
More information can be obtained from Kuiper, Hoijtink, and Silvapulle (2011) and Kuiper, Hoijtink, and Silvapulle (unpublished).
}
}
\value{an object of class orgls}
\references{
Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
\bold{98}, 495--501.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the Order-Restricted Information Criterion for Multivariate Normal Linear Models. \emph{Journal of Statistical Planning and Inference}, \bold{142}, 2454-2463. doi:10.1016/j.jspi.2012.03.007.

Kuiper R.M. and Hoijtink H. (submitted). A Fortran 90 Program for the Generalization of the Order-Restricted Information Criterion. Journal of Statictical Software.
}
\author{Daniel Gerhard and Rebecca M. Kuiper}


\seealso{\code{\link{solve.QP}}, \code{\link{goric}}}

\examples{
# generating example data
library(mvtnorm)
# group means
m <- c(0,5,5,7)
# compound symmetry structure of residuals
# (10 individuals per group, rho=0.7) 
cormat <- kronecker(diag(length(m)*10), matrix(0.7, nrow=length(m), ncol=length(m)))
diag(cormat) <- 1
# different variances per group
sds <- rep(c(1,2,0.5,1), times=10*length(m))
sigma <- crossprod(diag(sds), crossprod(cormat, diag(sds)))
response <- as.vector(rmvnorm(1, rep(m, times=10*length(m)), sigma=sigma))
dat <- data.frame(response,
  grp=rep(LETTERS[1:length(m)], times=10*length(m)),
  ID=as.factor(rep(1:(10*length(m)), each=length(m))))

## set of gls models:
# unconstrained model
m1 <- orgls(response ~ grp-1, data = dat,
  constr=rbind(c(0,0,0,0)), rhs=0, nec=0,
  weights=varIdent(form=~1|grp),
  correlation=corCompSymm(form=~1|ID))

# simple order
m2 <- orgls(response ~ grp-1, data = dat,
  constr=rbind(c(-1,1,0,0),c(0,-1,1,0),c(0,0,-1,1)), rhs=c(0,0,0), nec=0,
  weights=varIdent(form=~1|grp),
  correlation=corCompSymm(form=~1|ID))

# equality constraints
m3 <- orgls(response ~ grp-1, data = dat,
  constr=rbind(c(-1,1,0,0),c(0,-1,1,0),c(0,0,-1,1)), rhs=c(0,0,0), nec=3,
  weights=varIdent(form=~1|grp),
  correlation=corCompSymm(form=~1|ID))
}

\keyword{models}
