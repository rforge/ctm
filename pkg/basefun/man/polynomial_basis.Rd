\name{polynomial_basis}
\alias{polynomial_basis}
\title{
    Polynomial Basis Functions
}
\description{
    Basis functions defining a polynomial
}
\usage{
polynomial_basis(var, coef, ui = NULL, ci = NULL)
}
\arguments{
  \item{var}{a \code{\link[variables]{numeric_var}} object}
  \item{coef}{a logical defining the order of the polynomial}
  \item{ui}{a matrix defining constraints}
  \item{ci}{a vector defining constraints}
}
\details{
  \code{polynomial_basis} returns a function for the evaluation of
  the basis functions with corresponding \code{model.matrix} and \code{predict}
  methods.
}
\examples{

  ### set-up basis of order 3 ommiting the quadratic term
  pb <- polynomial_basis(numeric_var("x", support = c(0, pi)), 
                         coef = c(TRUE, TRUE, FALSE, TRUE))

  ### generate data + coefficients
  x <- as.data.frame(mkgrid(pb, n = 100))
  cf <- c(1, 2, 0, 1.75)

  ### evaluate basis (in two equivalent ways)
  pb(x[1:10,,drop = FALSE])
  model.matrix(pb, data = x[1:10, ,drop = FALSE])

  ### evaluate and plot polynomial defined by
  ### basis and coefficients
  plot(x$x, predict(pb, newdata = x, coef = cf), type = "l")

}
