\name{trafotree}
\alias{trafotree}
\title{Transformation Trees}
\description{
  Partitioned transformation models
}
\usage{
trafotree(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...)
}
\arguments{
  \item{object}{an object of class \code{\link[mlt]{ctm}} specifying the
                abstract model to be partitioned.}
  \item{parm}{parameters of \code{object} those corresponding score is
              used for finding partitions.}
  \item{mltargs}{arguments to \code{\link[mlt]{mlt}} for fitting the
                 transformation models.}
  \item{\dots}{arguments to \code{\link[partykit]{ctree}}, at least
               \code{formula} and \code{data}.}
}
\details{
  Conditional inference trees are used for partitioning likelihood-based transformation
  models as described in Hothorn and Zeileis (2017).
}
\value{
  An object of class \code{trafotree} with corresponding \code{plot}, \code{logLik} and
  \code{predict} methods.}
\references{
  Torsten Hothorn and Achim Zeileis (2017). Transformation Forests.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trees}
