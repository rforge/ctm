
library("tram")
library("mvtnorm")
set.seed(29)

ll <- numeric(50)

N <- 5000
p <- 3
X <- matrix(runif(N * p), ncol = p)
m1 <- 1 + X %*% c(2, 1, 1)
m2 <- 1 + X %*% c(1, 2, 1)
lb <- (off <- 0.5) + X %*% (cf <- c(0, 2, 0))
d <- data.frame(X)
Y <- matrix(NA, nrow = N, ncol = 2)
colnames(Y) <- c("Y1", "Y2")

cr <- numeric(N)
for (i in 1:N) {
    L <- diag(2)
    L[2,1] <- lb[i]
    Si <- solve(L) %*% t(solve(L))
#    Si <- diag(2)
#    Si[1,2] <- Si[2,1] <- .5
    cr[i] <- cov2cor(Si)[2,1]

    Y[i,] <- rmvnorm(1, mean = c(m1[i], m2[i]), sigma = Si)
}

d <- cbind(d, Y)
b1 <- as.mlt(Lm(Y1 ~ X1 + X2 + X3, data = d))
b2 <- as.mlt(Lm(Y2 ~ X1 + X2 + X3, data = d))

mm01 <- mmlt(b1, b2, formula = ~ X1 + X2 + X3, data = d, diag = FALSE)
mm02 <- mmlt(b2, b1, formula = ~ X1 + X2 + X3, data = d, diag = FALSE)

logLik(mm01) 
logLik(mm02)


x <- 0:4 / 4
nd <- expand.grid(X1 = x, X2 = x, X3 = x)
lhat <- off + as.matrix(nd) %*% cf
chat <- sapply(lhat, function(x) {
    L <- diag(2)
    L[2,1] <- x
    Si <- solve(L) %*% t(solve(L))
    cov2cor(Si)[2,1]
})
CR <- cbind(chat, coef(mm1, newdata = nd, type = "Cor"), 
            coef(mm2, newdata = nd, type = "Cor"))
pairs(CR)
abline(a = 0, b = 1)

predict(mm01, newdata = nd[1:5,], q = -2:2, 
  marginal = 1, type = "distribution")
predict(mm02, newdata = nd[1:5,], q = -2:2, 
  marginal = 2, type = "distribution")

predict(mm1, newdata = nd[1:5,], q = -2:2, 
  marginal = 1, type = "distribution")
predict(mm2, newdata = nd[1:5,], q = -2:2, 
  marginal = 2, type = "distribution")

c(coef(mm01, newdata = nd[1:5,], type = "Cor"))
c(coef(mm02, newdata = nd[1:5,], type = "Cor"))
c(coef(mm1, newdata = nd[1:5,], type = "Cor"))
c(coef(mm2, newdata = nd[1:5,], type = "Cor"))



## checking gradient for diag = TRUE
all.equal(c(numDeriv::grad(mm1$ll, mm2$par)),c(mm1$sc(mm2$par)), 
          check.attributes = FALSE, tol = 1e-4)
numDeriv::grad(mm1$ll, mm2$par)
mm1$sc(mm2$par)

all.equal(c(numDeriv::grad(mm2$ll, mm1$par)),c(mm2$sc(mm1$par)), 
          check.attributes = FALSE, tol = 1e-4)
numDeriv::grad(mm2$ll, mm1$par)
mm2$sc(mm1$par)

## evaluating logLik for a multiple of lambda
## only consider multiples of the off-diagonal elements?
logLik(mm1)
mm1$ll(mm1$par)
newpar <- c(mm1$pars$mpar, mm1$pars$cpar[, 1], 2*mm1$pars$cpar[, 2],
            mm1$pars$cpar[, 3])
mm1$ll(newpar)


logLik(mm01)
mm01$ll(mm01$par)
newpar <- c(mm01$pars$mpar, 2*mm01$pars$cpar)
mm01$ll(newpar)
