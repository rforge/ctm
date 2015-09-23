\name{factor_var}
\alias{factor_var}
\title{
    Unordered categorical variable
}
\description{
    Formal description of an unordered categorical variable
}
\usage{
factor_var(name, desc = NULL, levels)
}
\arguments{
  \item{name}{character, the name of the variable}
  \item{desc}{character, a description of what is measured}
  \item{levels}{character, the levels of the factor}
}
\details{
  A conceptual description of a (yet) unobserved unordered categorical variable.
}
\value{
  An object of class \code{factor\_var} inheriting from \code{var} with
  corresponding methods.
}
\examples{

  factor_var("eye", "eye color", c("blue", "brown", "green", "grey", "mixed"))

}
