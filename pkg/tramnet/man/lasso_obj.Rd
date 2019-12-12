
\name{lasso_obj}
\alias{lasso_obj}
\title{Lasso objective function for model based optimization}
\usage{
lasso_obj(object, minlambda = 0, maxlambda = 16, folds,
  noisy = FALSE, fold)
}
\arguments{
\item{object}{object of class \code{tramnet}}

\item{minlambda}{minimum value for lambda (default: \code{0})}

\item{maxlambda}{maximum value for lambda (default: \code{16})}

\item{folds}{self specified folds for cross validation (mainly for reproducibility
and comparability purposes)}

\item{noisy}{indicates whether folds for k-fold cross-validation should
be random for each iteration, leading to a noisy objective function
(default: \code{FALSE})}

\item{fold}{fold for cross validation}
}
\description{
Lasso objective function for model based optimization
}
