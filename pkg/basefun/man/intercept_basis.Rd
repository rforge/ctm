\name{intercept_basis}
\alias{intercept_basis}
\title{
    Intercept-Only Basis Function
}
\description{
    A simple intercept as basis function
}
\usage{
intercept_basis(ui = c("none", "increasing", "decreasing"), negative = FALSE)
}
\arguments{
  \item{ui}{a character describing possible constraints}
  \item{negative}{a logical indicating negative basis functions}
}
\details{
  \code{intercept_basis} returns a function for the evaluation of
  the basis functions with corresponding \code{model.matrix} and \code{predict}
  methods.
}
\examples{

  ### set-up basis
  ib <- intercept_basis()

  ### generate data + coefficients
  x <- as.data.frame(mkgrid(ib))

  ### 2 * 1 
  predict(ib, newdata = x, coef = 2)

}
