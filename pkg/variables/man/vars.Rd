\name{vars}
\alias{as.vars}
\alias{c.var}
\title{
    Multiple Abstract Descriptions
}
\description{
    Concatenate or generate multiple variable descriptions
}
\usage{
\method{c}{var}(...)
as.vars(object)
}
\arguments{
  \item{object}{an object}
  \item{\dots}{a list of variable objects}
}
\details{
  \code{c()} can be used to concatenate multiple variable objects; the corresponding
   generics also work for the resulting object. \code{as.vars()} tries to infer a 
  formal description from data.
}
\examples{

   f <- factor_var("x", levels = LETTERS[1:3])
   n <- numeric_var("y")

   fn <- c(f, n)
   variable.names(fn)
   support(fn)
   is.bounded(fn)
   mkgrid(fn, n = 9)

   as.vars(iris)
}
