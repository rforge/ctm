
library("mlt")

set.seed(28)
n <- 100
g <- gl(2, n)
y <- rnorm(length(g), mean = c(2, 1)[g], sd = c(.5, 1.5)[g])
mydata <- data.frame(g = g, y = y)

by <- polynomial_basis(c(1, 1), var = "y", support = range(y), ci = c(-Inf, 0))
bg1 <- as.basis(~ g - 1)

m1 <- model(by, inter = bg1, remove_intercept = TRUE)
fm1 <- mlt(m1, data = mydata)
logLik(fm1)
vcov(fm1)
cf1 <- coef(fm1)
1 / cf1[c(2, 4)] 
-cf1[c(1, 3)] / cf1[c(2,4)]

### no constraints here!
by2 <- polynomial_basis(c(1, 1), var = "y", support = range(y))
bg2 <- as.basis(~ g, remove_intercept = TRUE)
m2 <- model(by, inter = b(b1 = by2, b2 = bg2))
fm2 <- mlt(m2, data = mydata)
logLik(fm2)
vcov(fm2)
(cf2 <- coef(fm2))
c(cf1[1:2], cf1[3:4] - cf1[1:2])

1 / cf2[2]
1 / sum(cf2[c(2, 4)])

-cf2[1] / cf2[2]
-sum(cf2[c(1, 3)]) / sum(cf2[c(2, 4)])
