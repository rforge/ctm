\name{c.basis}
\alias{c.basis}
\title{
    Join Basis Functions
}
\description{
    Concatenate basis functions column-wise
}
\usage{
\method{c}{basis}(..., recursive = FALSE)
}
\arguments{
  \item{\dots}{named objects of class \code{basis}}
  \item{recursive}{always \code{FALSE}}
}
\details{
  \code{c()} joins the corresponding design matrices 
  column-wise, ie, the two functions defined by the two bases
  are added.
}
\examples{

  ### set-up Bernstein and log basis functions
  xv <- numeric_var("x", support = c(1, pi))
  bb <- Bernstein_basis(xv, order = 3, ui = "increasing")
  lb <- log_basis(xv, remove_intercept = TRUE)
  
  ### join them
  blb <- c(bern = bb, log = lb)

  ### generate data + coefficients
  x <- as.data.frame(mkgrid(blb, n = 100))
  cf <- c(1, 2, 2.5, 2.6, 2)

  ### evaluate bases
  model.matrix(blb, data = x[1:10, ,drop = FALSE])

  ### evaluate and plot function defined by
  ### bases and coefficients
  plot(x$x, predict(blb, newdata = x, coef = cf), type = "l")

  ### evaluate and plot first derivative of function
  ### defined by bases and coefficients
  plot(x$x, predict(blb, newdata = x, coef = cf, deriv = c(x = 1)), 
       type = "l")

}
