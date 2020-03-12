
\name{cvl_tramnet}
\alias{cvl_tramnet}
\title{Cross validation for \code{"tramnet"} models}
\usage{
cvl_tramnet(object, fold = 2, lambda = 0, alpha = 0, folds = NULL,
  fit_opt = FALSE)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{fold}{number of folds for cross validation}

\item{lambda}{values for lambda to iterate over}

\item{alpha}{values for alpha to iterate over}

\item{folds}{manually specify folds for comparison with other methods}

\item{fit_opt}{If \code{TRUE}, returns the full model evaluated at optimal hyper
parameters}
}
\value{
Returns out-of-sample logLik and coefficient estimates for
    corresponding folds and values of the hyperparameters as an object of
    class \code{"cvl_tramnet"}
}
\description{
k-fold cross validation for \code{"tramnet"} objects over a grid of
    the tuning parameters based on out-of-sample log-likelihood.
}
\examples{
\donttest{
set.seed(241068)
library(survival)
data("GBSG2", package = "TH.data")
X <- 1 * matrix(GBSG2$horTh == "yes", ncol = 1)
colnames(X) <- "horThyes"
GBSG2$surv <- with(GBSG2, Surv(time, cens))
m <- Coxph(surv ~ 1, data = GBSG2)
mt <- tramnet(model = m, x = X, lambda = 0, alpha = 0)
mc <- Coxph(surv ~ horTh, data = GBSG2)
cvl_tramnet(mt, fold = 2, lambda = c(0, 1), alpha = c(0, 1))
}
}
\author{
Lucas Kook
}
