\name{ordered_var}
\alias{ordered_var}
\title{
    Ordered Categorical Variable
}
\description{
    Formal description of an ordered categorical variable
}
\usage{
ordered_var(name, desc = NULL, levels)
}
\arguments{
  \item{name}{character, the name of the variable}
  \item{desc}{character, a description of what is measured}
  \item{levels}{character, the ordered levels of the factor}
}
\details{
  A conceptual description of a (yet) unobserved ordered categorical variable.
}
\value{
  An object of class \code{ordered\_var} inheriting from \code{var} with
  corresponding methods.
}
\examples{

  ordered_var("temp", "temperature", c("cold", "lukewarm", "warm", "hot"))

}
