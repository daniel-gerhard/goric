\name{constrSet}
\alias{constrSet}

\title{Generate Constraint Sets}
\description{Generate sets of constraint matrices (constr), right hand
  side elements, and numbers of equality constraints (nec) with a predefined structure}

\usage{
constrSet(n,
set=c("sequence","seqcontrol","lplateau","uplateau","downturn","williams"),
direction=c("increase","decrease"), base=1)
}
\arguments{
  \item{n}{a (possibly named) vector of sample sizes for each group.}
  \item{set}{character string defining the type of constraints; one of
    "sequence", "seqcontrol", "lplateau", "uplateau", or "downturn"}
  \item{direction}{direction of the inequality constraints, either
    "increase" or "decrease"}
  \item{base}{column of the constraint matrix representing a control
    group}
}
\value{a list with slots constr, rhs, and nec for each constraint definition}

\author{Daniel Gerhard and Rebecca M. Kuiper}


\seealso{\code{\link{orlm}}, \code{\link{constrSet}}}

\examples{
n <- c(10,20,30,40)
constrSet(n, set="sequence")
constrSet(n, set="seqcontrol")
constrSet(n, set="lplateau")
constrSet(n, set="uplateau")
constrSet(n, set="downturn")
constrSet(n, set="williams")
}

\keyword{misc}
