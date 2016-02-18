\name{access}
\alias{variable.names.var}
\alias{desc}
\alias{unit}
\alias{support}
\alias{bounds}
\alias{is.bounded}
\title{
    Accessor Functions
}
\description{
    Access properties of variable objects
}
\usage{
\method{variable.names}{var}(object, ...)
desc(object)
unit(object)
support(object)
bounds(object)
is.bounded(object)
}
\arguments{
  \item{object}{a variable object}
  \item{\dots}{additional arguments, currently not used}
}
\details{
  These generics have corresponding methods for \code{\link{factor_var}},
  \code{\link{ordered_var}} and \code{\link{numeric_var}} objects as well
  as for \code{vars} collections of those.
}
