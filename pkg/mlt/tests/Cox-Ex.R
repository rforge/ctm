
library("mlt")
library("survival")
set.seed(29)

### true dgp
rY <- function(n, ...) rexp(n, ...)
pY <- function(x, ...) pexp(x, ...)
dY <- function(x, ...) dexp(x, ...)

gf <- gl(3, 1)
g <- rep(gf, 100)
y <- rY(length(g), rate = (1:nlevels(g))[g])
mydata <- data.frame(y = y, g = g)

boxplot(y ~ g, data = mydata)

Bb <- Bernstein_basis(order = 5, support = c(0, max(y) + .1),
                      ui = "increasing", var = "y")
s <- as.basis(~ g, data = data.frame(g = gf), remove_intercept = TRUE)

m <- model(response = Bb, shifting = s, todist = "MinExtrVal")

(cf1 <- coef(opt <- mlt(m, data = mydata)))

library("lattice")
plot(opt, ~ y | g, what = "prob", plotfun = xyplot)

coef(cph <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata))

lb <- log_basis(varname = "y", support = c(0, max(y)))
m <- model(response = c(blog = lb, bBern = Bb), shifting = s, todist = "MinExtrVal")

(cf1 <- coef(opt <- mlt(m, data = mydata)))

nd <- samplefrom(opt, newdata = data.frame(g = gf), n = 100)

### check interval censoring!!!
(cf2 <- coef(opt2 <- mlt(m, data = nd)))

cf1 - cf2

yn <- generate(Bb, 50)$y
yn <- yn[yn > 0]

a <- predict(opt, newdata = data.frame(g = gf))

layout(matrix(1:4, ncol = 2))
plot(yn, a[[1]](yn, type = "trafo"), type = "l", col = "red")
lines(yn, log(yn))
plot(yn, 1 - a[[1]](yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[1])))
plot(yn, 1 - a[[2]](yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[2])))
plot(yn, 1 - a[[3]](yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[3])))

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE)), g = g)
coef(opt <- mlt(m, data = mydata))

coef(cph <- coxph(y ~ g, data = mydata))

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE), type = "left"), g = g)
coef(opt <- mlt(m, data = mydata))

Bb <- Bernstein_basis(order = 5, support = c(0, max(y + 1) + .1),
                      ui = "increasing", var = "y")

mydata <- data.frame(y = Surv(y, y + 1, sample(0:3, length(y), replace = TRUE), type = "interval"), 
                     g = g)
coef(opt <- mlt(m, data = mydata))

mydata <- data.frame(y = y, g = g)

Bb2 <- Bernstein_basis(order = 5, support = c(0, max(y) + .1), var = "y")
### basically works wo constraints (where data are); cf 2sample.R
#                       ui = "increasing", ci = -sqrt(.Machine$double.eps))

if (FALSE) {
m <- model(response = Bb, 
           interacting = b(b1 = Bb2, b2 = as.basis(~ g, remove_intercept = TRUE)),
           todist = "MinExtrVal")

coef(opt <- mlt(m, data = mydata))

coef(cph <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata))

a <- predict(opt, newdata = data.frame(g = gf))

layout(matrix(1:4, ncol = 2))
plot(yn, a[[1]](yn), type = "l", col = "red")
lines(yn, log(yn))
plot(yn, 1 - a[[1]](yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[1])))
plot(yn, 1 - a[[2]](yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[2])))
plot(yn, 1 - a[[3]](yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[3])))
}
