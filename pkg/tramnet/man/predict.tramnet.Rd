
\name{predict.tramnet}
\alias{predict.tramnet}
\title{predict method for class \code{"tramnet"}}
\usage{
\method{predict}{tramnet}(object, newdata = .get_tramnet_data(object),
  ...)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{newdata}{data used for prediction}

\item{...}{Additional arguments to \code{predict.ctm}}
}
\description{
predict method for class \code{"tramnet"}
}
\author{
Lucas Kook
}
\value{
Vector of predictions based on \code{object} evaluated at each row
of \code{newdata}
}