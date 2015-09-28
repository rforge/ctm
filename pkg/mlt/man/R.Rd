\name{R}
\alias{R}
\alias{R.Surv}
\alias{R.factor}
\alias{R.ordered}
\alias{R.numeric}
\title{
    Response variable
}
\description{
    Represent a possibly censored or truncated response variable
}
\usage{
R(object, ...)
\method{R}{numeric}(object = NA, cleft = NA, cright = NA,
                    tleft = NA, tright = NA, ...)
\method{R}{ordered}(object, cleft = NA, cright = NA, ...)
\method{R}{factor}(object, ...)
\method{R}{Surv}(object, ...)
}
\arguments{
  \item{object}{
  A vector of (conceptually) exact measurements}
  \item{cleft}{
  A vector of left borders of censored measurements}
  \item{cright}{
  A vector of right borders of censored measurements}
  \item{tleft}{
  A vector of left truncations}
  \item{tright}{
  A vector of right truncations}
  \item{\dots}{
  other arguments, ignored except for \code{tleft} and \code{tright} to
  \code{R.ordered}}
}
\details{
  \code{R} is basically an extention of \code{\link[survival]{Surv}} for the
  representation of arbitrarily censored or truncated measurements at any scale.
}
\examples{

 ### ordered factor
 R(gl(3, 3, labels = LETTERS[1:3]))

}
