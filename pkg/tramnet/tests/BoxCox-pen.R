# Tests for Box-Cox type regularized regression models

## Dependencies
library("tramnet")
options(digits = 3)

## Data
data("BostonHousing2", package = "mlbench")

## Set up and model fit
X <- model.matrix(cmedv ~ lstat + tax, data = BostonHousing2)[,-1]
m <- BoxCox(cmedv ~ 1, data = BostonHousing2)
mt <- tramnet(m, X, alpha = 0, lambda = 0)
m2 <- BoxCox(cmedv ~ lstat + tax, data = BostonHousing2)

max(abs(coef(mt, with_baseline = FALSE) -
        coef(m2, with_baseline = FALSE)))
logLik(mt)
logLik(m2)

## Test for additional inequality constraints on beta
mt <- tramnet(m, X, alpha = 0, lambda = 0, beta_const = matrix(0:1, nrow = 1))
m2 <- BoxCox(cmedv ~ lstat + tax, data = BostonHousing2, fixed = c("tax" = 0))

max(abs(coef(mt, with_baseline = FALSE) -
        coef(m2, with_baseline = FALSE)[-2]))
logLik(mt)
logLik(m2)
