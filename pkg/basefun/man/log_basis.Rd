\name{log_basis}
\alias{log_basis}
\title{
    Logarithmic Basis Function
}
\description{
    The logarithmic basis function
}
\usage{
log_basis(var, ui = c("none", "increasing", "decreasing"),
          remove_intercept = FALSE)
}
\arguments{
  \item{var}{a \code{\link[variables]{numeric_var}} object}
  \item{ui}{a character describing possible constraints}
  \item{remove_intercept}{a logical indicating if the intercept 
        term shall be removed}
}
\details{
  \code{log_basis} returns a function for the evaluation of
  the basis functions with corresponding \code{model.matrix} and \code{predict}
  methods.
}
\examples{

  ### set-up basis
  lb <- log_basis(numeric_var("x", support = c(0.1, pi)))

  ### generate data + coefficients
  x <- as.data.frame(mkgrid(lb, n = 100))

  ### 1 + 2 * log(x) 
  max(abs(predict(lb, newdata = x, coef = c(1, 2)) - (1 + 2 * log(x$x))))

}
