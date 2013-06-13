\name{constrMat}
\alias{constrMat}

\title{Generate Constraint Matrices}
\description{Generate a constraint matrix with a predefined structure}

\usage{
constrMat(n, type=c("monotone","control","average",
                    "laverage","uaverage","caverage"),
          base=1)
}
\arguments{
  \item{n}{a (possibly named) vector of sample sizes for each group}
  \item{type}{character string defining the type of constraints;
    one of
    "monotone", "control","average","laverage","uaverage", or "caverage"}
  \item{base}{column of the constraint matrix representing a control
    group (when type = "control")}
}
\value{a constraint matrix}

\author{Daniel Gerhard and Rebecca M. Kuiper}


\seealso{\code{\link{orlm}}, \code{\link{constrSet}}}

\examples{
n <- c(10,20,30,40)
constrMat(n, type="monotone")
constrMat(n, type="control", base=2)
constrMat(n, type="average")
constrMat(n, type="laverage")
constrMat(n, type="uaverage")
constrMat(n, type="caverage", base=2)
}

\keyword{misc}
