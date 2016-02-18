\name{Legendre_basis}
\alias{Legendre_basis}
\title{
    Legendre Basis Functions
}
\description{
    Basis functions defining a Legendre polynomial
}
\usage{
Legendre_basis(var, order = 2, ui = c("none", "increasing", "decreasing", 
                                      "cyclic", "positive", "negative"), ...)
}
\arguments{
  \item{var}{a \code{\link[variables]{numeric_var}} object}
  \item{order}{the order of the polynomial, one defines a linear function}
  \item{ui}{a character describing possible constraints}
  \item{\dots}{additional arguments passed to \code{\link[orthopolynom]{legendre.polynomials}}}
}
\details{
  \code{Legendre_basis} returns a function for the evaluation of
  the basis functions with corresponding \code{model.matrix} and \code{predict}
  methods.
}
\references{

  Rida T. Farouki (2012), The Bernstein Polynomial Basis: A Centennial Retrospective,
  \emph{Computer Aided Geometric Design}, \bold{29}(6), 379--419.
  http://dx.doi.org/10.1016/j.cagd.2012.03.001

}
\examples{

  ### set-up basis
  lb <- Legendre_basis(numeric_var("x", support = c(0, pi)), 
                       order = 3)

  ### generate data + coefficients
  x <- as.data.frame(mkgrid(lb, n = 100))
  cf <- c(1, 2, 2.5, 1.75)

  ### evaluate basis (in two equivalent ways)
  lb(x[1:10,,drop = FALSE])
  model.matrix(lb, data = x[1:10, ,drop = FALSE])

  ### evaluate and plot Legendre polynomial defined by
  ### basis and coefficients
  plot(x$x, predict(lb, newdata = x, coef = cf), type = "l")

}
