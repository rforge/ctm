\name{numeric_var}
\alias{numeric_var}
\title{
  Numeric variable
}
\description{
  Formal description of numeric variable
}
\usage{
numeric_var(name, desc = NULL, unit = NULL, support = c(0, 1), bounds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{character, the name of the variable}
  \item{desc}{character, a description of what is measured}
  \item{unit}{character, the measurement unit}
  \item{support}{the support of the measurements, see below}
  \item{bounds}{an interval defining the bounds of a real sample space}
}
\details{
  A numeric variable can be discrete (support is then the set of all
  possible values, either either integer or double; integers of length 2
  are interpreted as all integers inbetween) or continuous (support is
  a double of length 2 giving the support of the data).

  If a continuous variable is bounded, \code{bounds} defines the 
  corresponding interval.
}
\value{
  An object of class \code{numeric\_var} inheriting from \code{var} with
  corresponding methods.
}
\examples{

  numeric_var("age", "age of patient", "years", support = 25:75)

  numeric_var("time", "survival time", "days", support = 0:365)

  numeric_var("time", "survival time", "days", support = c(0.0, 365), bounds = c(0, Inf))

}
