### test for mcotram
library("cotram")
library("mvtnorm")
library("Matrix")
set.seed(123)
N <- 1000
S <- diag(2)
S[1, 2] <- S[2, 1] <- .9

x <- matrix(runif(N * 2), ncol = 2)

beta <- c(1, -1, -.5, .5)
beta <- rep(0, 4)

ly <- x %*% matrix(beta, nrow = 2) + rmvnorm(N, sigma = S)

my <- ceiling(max(exp(ly)))
ey <- exp(ly)
y <- ey

y[,1] <- unclass(oy1 <- cut(ey[,1], breaks = 0:my, ordered = TRUE))
y[,2] <- unclass(oy2 <- cut(ey[,2], breaks = 0:my, ordered = TRUE))

d <- data.frame(y = y, x = x, ey = ey, oy1 = oy1, oy2 = oy2)

b1 <- BoxCox(ey.1 ~ x.1 + x.2, data = d)
coef(b1)
logLik(b1)
b2 <- BoxCox(ey.2 ~ x.1 + x.2, data = d)
coef(b2)
logLik(b2)
p1 <- Polr(oy1 ~ x.1 + x.2, data = d, method = "probit")
coef(p1)
logLik(p1)
p2 <- Polr(oy2 ~ x.1 + x.2, data = d, method = "probit")
coef(p2)
logLik(p2)

u1 <- cotram(y.1 ~ 1, data = d, method = "probit")
u2 <- cotram(y.2 ~ 1, data = d, method = "probit")
uc1 <- mcotram(u1, u2, data = d)
uc2 <- mcotram(u2, u1, data = d)
logLik(uc1)
logLik(uc2)
uc1 <- mcotram(u1, u2, data = d, gr = FALSE)
uc2 <- mcotram(u2, u1, data = d, gr = FALSE)
logLik(uc1)
logLik(uc2)


u1 <- cotram(y.1 ~ 1, data = d, method = "probit", log_first = FALSE)
u2 <- cotram(y.2 ~ 1, data = d, method = "probit", log_first = FALSE)
uc1 <- mcotram(u1, u2, data = d)
uc2 <- mcotram(u2, u1, data = d)
logLik(uc1)
logLik(uc2)
uc1 <- mcotram(u1, u2, data = d, gr = FALSE)
uc2 <- mcotram(u2, u1, data = d, gr = FALSE)
logLik(uc1)
logLik(uc2)



m01 <- cotram(y.1 ~ x.1 + x.2, data = d, method = "probit", log_first =
FALSE)
coef(m01)
logLik(m01)
m02 <- cotram(y.2 ~ x.1 + x.2, data = d, method = "probit", log_first =
FALSE)
coef(m02)
logLik(m02)

## simple formula
mc01 <- mcotram(m01, m02, formula = ~ 1, data = d)
mc02 <- mcotram(m02, m01, formula = ~ 1, data = d)

coef(mc01)
coef(mc02)

logLik(mc01)
logLik(mc02)

mb1 <- mmlt(b1, b2, formula = ~ 1, data = d)
mb2 <- mmlt(b2, b1, formula = ~ 1, data = d)

coef(mb1)
coef(mb2)

logLik(mb1)
logLik(mb2)



## complex formula
mc1 <- mcotram(m1, m2, formula = ~ x.1 + x.2, data = d, diag = FALSE)
mc2 <- mcotram(m2, m1, formula = ~ x.1 + x.2, data = d, diag = FALSE)
logLik(mc1)
logLik(mc2)

mc1d <- mcotram(m1, m2, formula = ~ x.1 + x.2, data = d, diag = TRUE)
mc2d <- mcotram(m2, m1, formula = ~ x.1 + x.2, data = d, diag = TRUE)
logLik(mc1d)
logLik(mc2d)

## checking gradient for diag = TRUE
all.equal(c(numDeriv::grad(mc1d$ll, mc2d$par)),c(mc1d$sc(mc2d$par)), 
          check.attributes = FALSE, tol = 1e-6)
numDeriv::grad(mc1d$ll, mc2d$par)
mc1d$sc(mc2d$par)

all.equal(c(numDeriv::grad(mc2d$ll, mc1d$par)),c(mc2d$sc(mc1d$par)), 
          check.attributes = FALSE, tol = 1e-6)
numDeriv::grad(mc2d$ll, mc1d$par)
mc2d$sc(mc1d$par)

## now with continuous mmlt
# library("tram")
# m1 <- as.mlt(BoxCox(y.1 ~ x.1 + x.2, data = d))
# m2 <- as.mlt(BoxCox(y.2 ~ x.1 + x.2, data = d))
# 
# mc1 <- mmlt(m1, m2, formula = ~ 1, data = d)
# mc2 <- mmlt(m2, m1, formula = ~ 1, data = d)
# 
# coef(m1)
# coef(mc1)
# coef(m2)
# coef(mc2)
# 
# logLik(mc1)
# logLik(mc2)
# 
# ### different formula in multivariate model:
# mc1 <- mmlt(m1, m2, formula = ~ x.1 + x.2, data = d)
# mc2 <- mmlt(m2, m1, formula = ~ x.1 + x.2, data = d)
# 
# coef(m1)
# coef(mc1)
# coef(m2)
# coef(mc2)
# 
# logLik(mc1)
# logLik(mc2)


## testing predict performance
predict(m1, newdata = d[1:3, ], type = "distribution")
predict(mc1, marginal = 1, newdata = d[1:3, ], type = "distribution")
predict(mc2, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(m2, newdata = d[1:3, ], type = "distribution")
predict(mc1, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(mc2, marginal = 1, newdata = d[1:3, ], type = "distribution")

head(coef(mc1, type = "Lambda"))
head(coef(mc2, type = "Lambda"))
head(coef(mc1, type = "Corr"))
head(coef(mc2, type = "Corr"))



predict(m1, newdata = d[1:3, ], type = "distribution")
predict(mc1d, marginal = 1, newdata = d[1:3, ], type = "distribution")
predict(mc2d, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(m2, newdata = d[1:3, ], type = "distribution")
predict(mc1d, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(mc2d, marginal = 1, newdata = d[1:3, ], type = "distribution")

head(coef(mc1d, type = "Corr"))
head(coef(mc2d, type = "Corr"))


################ dimension = 3 ################ 
set.seed(123)
N <- 1000
S <- diag(3)
S[lower.tri(S)] <- S[upper.tri(S)] <- .5

x <- matrix(runif(N*2), ncol = 2)
m <- x %*% c(-1, 1)

y <- x %*% matrix(c(1, -1, -.5, .5, -.2, .2), nrow = 2) + rmvnorm(N, sigma = S)

y[,1] <- unclass(cut(exp(y[,1]), breaks = c(0:10, Inf)))
y[,2] <- unclass(cut(exp(y[,2]), breaks = c(0:10, Inf)))
y[,3] <- unclass(cut(exp(y[,3]), breaks = c(0:10, Inf)))


d <- data.frame(y = y, x = x)

m1 <- cotram(y.1 ~ x.1 + x.2, data = d, method = "probit")
m2 <- cotram(y.2 ~ x.1 + x.2, data = d, method = "probit")
m3 <- cotram(y.3 ~ x.1 + x.2, data = d, method = "probit")

## simple formula
mc01 <- mcotram(m1, m2, m3, formula = ~ 1, data = d)
mc02 <- mcotram(m2, m3, m1, formula = ~ 1, data = d)

m1$coef
coef(mc01)
m2$coef
coef(mc02)

logLik(mc01)
logLik(mc02)

## complex formula
mc1 <- mcotram(m1, m2, m3, formula = ~ x.1 + x.2, data = d, diag = FALSE)
mc2 <- mcotram(m2, m3, m1, formula = ~ x.1 + x.2, data = d, diag = FALSE)
logLik(mc1)
logLik(mc2)

library(Matrix)
mc1d <- mcotram(m1, m2, m3, formula = ~ x.1 + x.2, data = d, diag = TRUE)
mc2d <- mcotram(m2, m3, m1, formula = ~ x.1 + x.2, data = d, diag = TRUE)
logLik(mc1d)
logLik(mc2d)

## checking gradient for diag = TRUE
all.equal(c(numDeriv::grad(mc1d$ll, mc2d$par)),c(mc1d$sc(mc2d$par)), 
          check.attributes = FALSE, tol = 1e-6)
numDeriv::grad(mc1d$ll, mc2d$par)
mc1d$sc(mc2d$par)

all.equal(c(numDeriv::grad(mc2d$ll, mc1d$par)),c(mc2d$sc(mc1d$par)), 
          check.attributes = FALSE, tol = 1e-6)
numDeriv::grad(mc2d$ll, mc1d$par)
mc2d$sc(mc1d$par)

## now with continuous mmlt
# library("tram")
# m1 <- as.mlt(Lm(y.1 ~ x.1 + x.2, data = d))
# m2 <- as.mlt(Lm(y.2 ~ x.1 + x.2, data = d))
# 
# mc1 <- mmlt(m1, m2, formula = ~ 1, data = d)
# mc2 <- mmlt(m2, m1, formula = ~ 1, data = d)
# 
# coef(m1)
# coef(mc1)
# coef(m2)
# coef(mc2)
# 
# logLik(mc1)
# logLik(mc2)
# 
# ### different formula in multivariate model:
# mc1 <- mmlt(m1, m2, formula = ~ x.1 + x.2, data = d)
# mc2 <- mmlt(m2, m1, formula = ~ x.1 + x.2, data = d)
# 
# coef(m1)
# coef(mc1)
# coef(m2)
# coef(mc2)
# 
# logLik(mc1)
# logLik(mc2)


## testing predict performance
predict(m1, newdata = d[1:3, ], type = "distribution")
predict(mc1, marginal = 1, newdata = d[1:3, ], type = "distribution")
predict(mc2, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(m2, newdata = d[1:3, ], type = "distribution")
predict(mc1, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(mc2, marginal = 1, newdata = d[1:3, ], type = "distribution")

head(coef(mc1, type = "Lambda"))
head(coef(mc2, type = "Lambda"))
head(coef(mc1, type = "Corr"))
head(coef(mc2, type = "Corr"))



predict(m1, newdata = d[1:3, ], type = "distribution")
predict(mc1d, marginal = 1, newdata = d[1:3, ], type = "distribution")
predict(mc2d, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(m2, newdata = d[1:3, ], type = "distribution")
predict(mc1d, marginal = 2, newdata = d[1:3, ], type = "distribution")
predict(mc2d, marginal = 1, newdata = d[1:3, ], type = "distribution")

head(coef(mc1d, type = "Corr"))
head(coef(mc2d, type = "Corr"))

