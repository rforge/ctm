
library("mlt")
library("truncreg")

set.seed(29)
n <- 1000	
x <- runif(n, max = 2 * pi)
y <- rnorm(n, 2*x + 1, .25)
d <- data.frame(y = y, x = x)
## truncated response
tr <- 2.5
d$yt <- ifelse(d$y > tr, d$y, NA)
d$trunc_left <- tr

tmod <- truncreg(yt ~ x, data = d, point = tr, direction = "left")
coef(tmod)
logLik(tmod)
vcov(tmod)

## MLT
fm <- as.formula("~ x")
yb <- polynomial_basis(c(TRUE, TRUE), var = "yt", ci = c(-Inf, 0))
m <- model(yb, shifting = fm, todistr = "Normal")
mltmod <- mlt(m, data = d[!is.na(d$yt),,drop = FALSE], trunc = c("left" = "trunc_left"))
(cf <- coef(mltmod))
c(-cf[1] / cf[2], -cf[3] / cf[2], 1 / cf[2])
logLik(mltmod)
vcov(mltmod)

library("numDeriv")

solve(hessian(mltmod$loglik, coef(mltmod)))