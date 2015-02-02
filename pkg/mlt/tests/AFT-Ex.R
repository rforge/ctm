
library("mlt")
set.seed(29)

### true dgp
rY <- function(n, ...) rexp(n, ...)
pY <- function(x, ...) pexp(x, ...)
dY <- function(x, ...) dexp(x, ...)

gf <- gl(3, 1)
g <- rep(gf, 10000)
y <- rY(length(g), rate = (1:nlevels(g))[g])
mydata <- data.frame(y = y, g = g)

boxplot(y ~ g, data = mydata)

Bb <- log_basis(var = "y", support = range(y))

Bx <- as.basis(~ g, remove_intercept = FALSE)
m <- model(Bb, shifting = Bx)

coef(opt <- mlt(m, data = mydata, fixed = c("log(y)" = 1),
                todist = "MinExtrVal"))

coef(aft <- survreg(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata,
                    dist = "exponential"))

yn <- generate(Bb, 50)$y

a1 <- predict(opt, newdata = data.frame(g = gf[1]), type = "trafo")
a2 <- predict(opt, newdata = data.frame(g = gf[2]), type = "trafo")
a3 <- predict(opt, newdata = data.frame(g = gf[3]), type = "trafo")

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE)), g = g)
coef(opt <- mlt(m, data = mydata, fixed = c("log(y)" = 1),
                todist = "MinExtrVal"))

coef(aft <- survreg(y ~ g, data = mydata, dist = "expo"))

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE), type = "left"), g = g)
coef(opt <- mlt(m, data = mydata,
                todist = "MinExtrVal"))

mydata <- data.frame(y = Surv(y, y + 1, sample(0:3, length(y), replace = TRUE), type = "interval"), 
                     g = g)
coef(aft <- mlt(m, data = mydata,
                todist = "MinExtrVal"))

