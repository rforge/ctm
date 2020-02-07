\name{score_test}
\alias{score_test}
\alias{score_test.tram}
\title{
  Transformation Score Tests and Confidence Intervals
}
\description{
  P-values and confidence intervals for parameters in linear transformation
  models obtained from by the score test principle
}
\usage{
score_test(object, ...)
\method{score_test}{tram}(object, parm = names(coef(object)), 
    alternative = c("two.sided", "less", "greater"), nullvalue = 0, 
    confint = TRUE, level = .95, Taylor = FALSE, maxsteps = 25, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{tram}}}
  \item{parm}{a vector of names of parameters to be tested.
              These parameters must be present in \code{object}.}
  \item{alternative}{a character string specifying the alternative hypothesis,
          must be one of \code{"two.sided"} (default), \code{"greater"} 
          or \code{"less"}.}
  \item{nullvalue}{a number specifying an optional parameter used to form the
          null hypothesis.}
  \item{confint}{a logical indicating whether a confidence interval should be
          computed. Score confidence intervals are computed by default. A
          1st order Taylor approximation to the Score statistc is used with
          \code{Taylor = TRUE} (in case numerical inversion of the score
          statistic fails, Wald confidence intervals relying from this approximation are
          returned).}
  \item{level}{the confidence level.}
  \item{Taylor}{a logical requesting the use of a 1st order Taylor
                approximation when inverting the score statistic.}
  \item{maxsteps}{number of function evaluations when inverting the score statistic
                  for computing confidence intervals.}
  \item{\dots}{additional arguments, currently ignored.}
}
\details{
  Score tests and confidence intervals for the parameters in the linear
  predictor of \code{object} are computed. These parameters must be present
  in \code{object}.
}
\value{
  An object of class \code{htest} or a list thereof. See \code{\link{Coxph}}
  for an example. A corresponding permutation test for parameters in a
  transformation models is available in
  \code{\link[tram]{perm_test}}.
}
