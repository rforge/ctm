# Tests for Lm models

## Dependencies
library("tramnet")
options(digits = 3)

## Data
dat <- data.frame(y = runif(100), s = factor(rep(c(1, 2), each = 50)))
x <- matrix(rnorm(100 * 20, mean = 0, sd = 1), nrow = 100)
colnames(x) <- paste0("X", 1:20)
y2 <- Lm(y | 0 + s ~ 1, data = dat)
mod2 <- tramnet(y2, x, lambda = 8, alpha = 1)

## --- Truncated, stratified
dat$ytr <- R(dat$y, tleft = c(-Inf, -2, rep(0, length(dat$y)-2)))
y3 <- Colr(ytr | 0 + s ~ 1, data = dat, order = 4)
## TODO: Implement truncation

