
library("tram")
library("sandwich")
library("survival")

### scores for fixed parameters
cars$int <- 1
a0 <- Lm(dist ~ speed, data = cars)
a1 <- Lm(dist ~ speed + int, data = cars, fixed = c("int" = 0))

s0 <- estfun(a0)
s1 <- estfun(a1)
s2 <- estfun(a1, parm = coef(as.mlt(a1), fixed = TRUE))

stopifnot(all.equal(s0[,"(Intercept)"], -s2[,"int"]))
stopifnot(all.equal(s1[,"(Intercept)"], -s2[,"int"]))
stopifnot(all.equal(s2[,"(Intercept)"], -s2[,"int"]))

### log_first for count data (by Sandra Siegfried)
### use: y + 1; log_first = TRUE, support = c(1, ...), bounds = c(1, ...)
library("MASS")
set.seed(29)
Nsim <- 100
b1 <- - 4.5
b0 <- 5
theta <- 2
h <- qlogis(ppois(0:100, lambda = 5))
dgp <- function(n = 4000){
  x <- runif(n, min = 0, max = 1) 
  log.mu <- b0 + b1 * x
  h.m <- matrix(h, nrow = length(h), ncol = length(x)) 
  p <- (plogis(t(h.m) - b1 * x) - runif(n))^2
  y <- max.col(-p) - 1
  m <- Colr(y ~ x, data = data.frame(y = y, x = x), 
            bounds = c(0L, Inf), fixed = c("x" = -b1),
            support = c(0L, floor(quantile(y, .9))), order = 10)
  
  ret <- data.frame(x = x, y = as.integer(y))
  attr(ret, "mC") <- m
  ret
}
d <- dgp()
d$y.p1 <- d$y + 1L

m1 <- Colr(y ~ x, data = d,
                  support = c(0L, as.numeric(max(d$y))),
                  bounds = c(0L, Inf),
                  order = 10)

try(m2a <- Colr(y ~ x, data = d,
                  support = c(1L, as.numeric(max(d$y))),
                  bounds = c(1L, Inf),
                  log_first = TRUE,
                  order = 10))

try(m2b <- Colr(y ~ x, data = d,
                  support = c(1L, as.numeric(max(d$y))),
                  bounds = c(0L, Inf),
                  log_first = TRUE,
                  order = 10))

m3 <- Colr(y.p1 ~ x, data = d,
                  support = c(1L, as.numeric(max(d$y.p1))),
                  bounds = c(1L, Inf),
                  log_first = TRUE,
                  order = 10)

l1 <- m1$logliki(coef(as.mlt(m1)), rep(1, nrow(d)))
l3 <- m3$logliki(coef(as.mlt(m3)), rep(1, nrow(d)))
stopifnot(cor(l1, l3) > .96)

### normal CDF approximation by Matic et al (2018)
x <- -600:600 / 100
p1 <- .Call("R_pnormMRS", x)
p2 <- pnorm(x)
stopifnot(max(abs(p1 - p2)) < 5.79e-6) ### max absolute error

### 0.3-0
data("GBSG2", package = "TH.data")
### this model included an additional intercept
m1 <- Coxph(Surv(time, cens) | menostat:tgrade ~ horTh, data = GBSG2)
m2 <- Coxph(Surv(time, cens) | 0 + menostat:tgrade ~ horTh, data = GBSG2)
stopifnot(max(abs(coef(as.mlt(m1)) - coef(as.mlt(m2)))) < 
          sqrt(.Machine$double.eps))
