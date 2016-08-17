
library("mlt")
set.seed(29)

n <- 2500
y <- rchisq(n, df = 1)

vy <- numeric_var("y", support = quantile(y, prob = c(.1, .9)), bounds = c(0, 10))
by <- Bernstein_basis(vy, order = 10, ui = "increasing")

m1 <- mlt(ctm(by), data = data.frame(y = y))

if (require("BDR")) {

mydf <- BDR(data.frame(y = y), as.int = "y", nmax = 50, total = TRUE)

mdf <- attr(mydf, "levels")

m2 <- mlt(ctm(by), data = mdf, weights = weights(mdf))

cf1 <- coef(m1)
cf2 <- coef(m2)

print(range(cf1/cf2))

}