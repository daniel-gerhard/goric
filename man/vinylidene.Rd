\name{vinylidene}
\alias{vinylidene}
\docType{data}
\title{Effect of vinylidene fluoride on liver cancer}
\description{Real data which are available on page 10 of Silvapulle and
  Sen (2005) and in a report prepared by Litton Bionetics Inc in 1984.
  
  These data were used in an experiment to find out whether
  vinylidene fluoride gives rise to liver damage. Since increased levels
  of serum enzyme are inherent in liver damage, the focus is on whether
  enzyme levels are affected by vinylidene fluoride.
  The variable of interest is the serum enzyme level. Three types of
  enzymes are inspected, namely SDH, SGOT, and SGPT. To study whether
  vinylidene fluoride has an influence on the three serum enzymes, four
  dosages of this substance are examined. In each of these four treatment groups, ten male Fischer-344 rats received the substance. 
}
\usage{data(vinylidene)}
\format{
  A data frame with 40 observations on the following 4 variables.
  \describe{
    \item{\code{SDH}}{serum enzyme level of enzyme type SDH.}
    \item{\code{SGOT}}{serum enzyme level of enzyme type SGOT.}
    \item{\code{SGPT}}{serum enzyme level of enzyme type SGPT.}
    \item{\code{dose}}{factor with 4 levels (d1-d4) representing the 4
      vinylidene fluoride concentrations.}
  }
}
\references{Silvapulle MJ, Sen PK (2005). \emph{Constrained Statistical Inference}. New Jersey: Wiley.}
\keyword{datasets}
