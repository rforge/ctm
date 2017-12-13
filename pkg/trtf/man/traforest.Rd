\name{traforest}
\alias{traforest}
\title{Transformation Forests}
\description{
  Partitioned and aggregated transformation models
}
\usage{
traforest(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...)
}
\arguments{
  \item{object}{an object of class \code{\link[mlt]{ctm}} specifying the
                abstract model to be partitioned.}
  \item{parm}{parameters of \code{object} those corresponding score is
              used for finding partitions.}
  \item{mltargs}{arguments to \code{\link[mlt]{mlt}} for fitting the
                 transformation models.}
  \item{\dots}{arguments to \code{\link[partykit]{cforest}}, at least
               \code{formula} and \code{data}.}
}
\details{
  Conditional inference trees are used for partitioning likelihood-based transformation
  models as described in Hothorn and Zeileis (2017). The method can be seen
  in action in Hothorn (2018) and the corresponding code is available as
  \code{demo("BMI")}. 
}
\value{
  An object of class \code{traforest} with corresponding \code{logLik} and
  \code{predict} methods.
}
\references{
  Torsten Hothorn and Achim Zeileis (2017). Transformation Forests.
  \url{https://arxiv.org/abs/1701.02110}.

  Torsten Hothorn (2018). Top-Down Transformation Choice. \emph{Statistical
  Modelling}, \url{https://arxiv.org/abs/1706.08269}.
}
\keyword{trees}
