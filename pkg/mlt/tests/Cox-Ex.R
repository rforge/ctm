
library("mlt")
set.seed(29)

### true dgp
rY <- function(n, ...) rexp(n, ...)
pY <- function(x, ...) pexp(x, ...)
dY <- function(x, ...) dexp(x, ...)

gf <- gl(3, 1)
g <- rep(gf, 30)
y <- rY(length(g), rate = (1:nlevels(g))[g])
mydata <- data.frame(y = y, g = g)

boxplot(y ~ g, data = mydata)

Bb <- Bernstein_basis(order = 10, support = c(0, max(y) + .1),
                      constraint = "increasing", var = "y")

opt <- mlt(response = "y", 
           data = mydata,
           bresponse = Bb,
           bshift = ~ g,
           dist = Distr("MinExtrVal"))

opt$par

coef(cph <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata))

yn <- generate(Bb, 50)

a1 <- predict(opt, newdata = data.frame(g = gf[1]), type = "trafo")
a2 <- predict(opt, newdata = data.frame(g = gf[2]), type = "trafo")
a3 <- predict(opt, newdata = data.frame(g = gf[3]), type = "trafo")

layout(matrix(1:4, ncol = 2))
plot(yn, a1(yn), type = "l", col = "red")
lines(yn, log(yn))
plot(yn, 1 - a1(yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[1])))
plot(yn, 1 - a2(yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[2])))
plot(yn, 1 - a3(yn, type = "prob"), type = "l", col = "red", ylim = c(0, 1))
lines(survfit(cph, newdata = data.frame(g = gf[3])))
