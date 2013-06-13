\name{goric-package}
\alias{goric-package}
\docType{package}
\title{
Generalized Order-Restricted Information Criterion (GORIC) value 
for a set of hypotheses in multivariate regression models
}
\description{
Computes the Generalized Order-Restricted Information Criterion (GORIC) values 
for a set of hypotheses in multivariate regression models. The GORIC is a modification of the AIC (Akaike, 1973) and the ORIC (Anraku, 1999) such that it can be applied to a general form of order-restricted hypotheses in multivariate normal linear models (Kuiper, Hoijtink, and Silvapulle, 2011; Kuiper, Hoijtink, and Silvapulle, unpublished).
}
\details{
\tabular{ll}{
Package: \tab goric\cr
Type: \tab Package\cr
Version: \tab 0.0-8\cr
Date: \tab 2013-03-26\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

In a \eqn{t}-variate regression with \eqn{k} predictors (including an optional intercept),
the order-restricted hypotheses should be of the form: 
\eqn{Constr \theta \geq rhs}, 
with 
\eqn{\theta} a vector of length \eqn{tk}, where the first \eqn{k} elements belong to first 
dependent variable, ..., and the last \eqn{k} elements to the last dependent variable;
\eqn{rhs} is a vector of length \eqn{c}; and 
\eqn{Constr} a \eqn{c \times tk} matrix of full rank.

Thera are two requirements:
\enumerate{
\item The first \eqn{nec} constraints must be the equality contraints 
(i.e., \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}) 
and the remaing ones the inequality contraints
(i.e., \eqn{Constr[nec+1:c, 1:tk] \theta \geq rhs[nec+1:c]}).
\item When \eqn{rhs} is not zero, \eqn{Constr} should be of full rank 
(after discarding redundant restrictions).
}
More information can be obtained from Kuiper, Hoijtink, and Silvapulle (2011), Kuiper, Hoijtink, and Silvapulle (2012), and Kuiper and Hoijtink (unpublished).


\code{orlm} renders the order-restricted maximum likelihood estimates
(i.e., the maximum likelihood estimates subject to the hypothesis of interest) 
and the corresponding log likelihood for the hypothesis of interest 
(defined by \eqn{Constr}, \eqn{rhs}, and \eqn{nec}).
Additionally it gives the (unconstrained) maximum likelihood estimates and
the active contraints.

\code{goric} gives the order-restricted log likelihood, the penalty of
the GORIC, the GORIC values, and the GORIC weights for a set of
hypotheses (orlm objects), where the penalty is based on \eqn{iter} iterations.

The hypothesis with the lowest GORIC value is the preferred one.

The GORIC weights reflect the support of each hypothesis in the set.
To compare two hypotheses (and not one to the whole set), one should examine 
the ratio of the two corresponding GORIC weights.

To safequard for weak hypotheses (i.e., hypotheses not supported by teh data), 
one should include a model with no constraints (the so-called unconstrained model).
More information can be obtained from Kuiper, Hoijtink, and Silvapulle (2011) and 
Kuiper, Hoijtink, and Silvapulle (unpublished).
}
\author{
Daniel Gerhard and Rebecca M. Kuiper

Maintainer: Daniel Gerhard <gerhard@biostat.uni-hannover.de>
}
\references{
Akaike, H. (1973). Information theory and an extension of the maximum likelihood principle. \emph{In Proc. 2nd Int. Symp. Information Theory}, Ed. B. N. Petrov 
and F. Csaki, pp. 267/226-81. Budapest: Akademiai kiado.

Anraku K. (1999). An Information Criterion for Parameters under a Simple Order
Restriction. \emph{Biometrika}, \bold{86}, 141--152.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
\bold{98}, 495--501.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the 
Order-Restricted Information Criterion for Multivariate Normal Linear Models. \emph{Journal of Statistical Planning and Inference}, 142, 2454-2463. doi:10.1016/j.jspi.2012.03.007.

Kuiper R.M. and Hoijtink H. (submitted to Journal of Statictical Software).
A Fortran 90 Program for the Generalization of the Order-Restricted Information Criterion.
}

\keyword{package}

