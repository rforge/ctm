
\name{mbo_tramnet}
\alias{mbo_tramnet}
\title{Model based optimization for regularized transformation models}
\usage{
mbo_tramnet(object, fold = 2, n_design = 5, n_iter = 5,
  minlambda = 0, maxlambda = 16, minalpha = 0, maxalpha = 1,
  folds = NULL, learner = "regr.km", pred.type = "se",
  opt_crit = makeMBOInfillCritEI(), noisy = FALSE,
  obj_type = c("lasso", "ridge", "elnet"), verbose = TRUE, ...)
}
\arguments{
\item{object}{object of class \code{tramnet}}

\item{fold}{fold for cross validation}

\item{n_design}{results in \code{n_design} times the number of tuning parameters
rows for the initial design matrix based on a random latin hypercube design}

\item{n_iter}{number of iterations in the model based optimization procedure}

\item{minlambda}{minimum value for lambda (default: \code{0})}

\item{maxlambda}{maximum value for lambda (default: \code{16})}

\item{minalpha}{minimum value for alpha (default: \code{0})}

\item{maxalpha}{maximum value for alpha (default: \code{1})}

\item{folds}{self specified folds for cross validation (mainly for reproducibility
and comparability purposes)}

\item{learner}{type of leaner used for the optimization (default: \code{"regr.km"})}

\item{pred.type}{prediction type of the learner (default: \code{"se"})}

\item{opt_crit}{optimization criterion, default: expected improvement}

\item{noisy}{indicates whether folds for k-fold cross-validation should
be random for each iteration, leading to a noisy objective function
(default: \code{FALSE})}

\item{obj_type}{objective type, one of \code{"lasso"}, \code{"ridge"} or \code{"elnet"}}

\item{verbose}{toggle for a verbose output (default: \code{TRUE})}

\item{...}{additional arguments are ignored}
}
\description{
Uses model based optimization to find the optimal tuning parameter(s)
    in a regularized transformation model based on cross-validated log-likelihoods.
    Here the tramnet package makes use of the mlrMBO interface for Bayesian
    Optimization in machine learning problems to maximize the cv-logLik
    as a black-box function of the tuning parameters alpha and lambda.
}
\examples{
\donttest{
set.seed(24101968)
data("Prostate", package = "lasso2")
x <- scale(model.matrix(lpsa ~ . - 1, data = Prostate))
y <- Prostate$lpsa
m0 <- BoxCox(lpsa ~ 1, data = Prostate)
mt <- tramnet(m0, x = x, lambda = 0, alpha = 1, check_dcp = FALSE)
(run_lasso <- mbo_tramnet(object = mt, n_design = 2, n_iter = 2,
                          obj_type = "lasso"))
}
}
\value{
returns an object of class \code{"MBOSingleObjResult"} which is
documented in \code{\link[mlrMBO]{mbo}}
}
