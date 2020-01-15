# Tests for Box-Cox type regularized regression models

## Dependencies
## IGNORE_RDIFF_BEGIN
library("tramnet")
## IGNORE_RDIFF_END
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
m2 <- BoxCox(cmedv ~ lstat + tax, data = BostonHousing2, constraints = c("tax >= 0"))
lhs <- attr(model.matrix(m2), "constraint")$ui
rhs <- attr(model.matrix(m2), "constraint")$ci
mt <- tramnet(m, X, alpha = 0, lambda = 0, 
              constraints = list(lhs, rhs))

max(abs(coef(mt, with_baseline = FALSE) -
        coef(m2, with_baseline = FALSE)[-2])) < 1e-5
logLik(mt)
logLik(m2)
