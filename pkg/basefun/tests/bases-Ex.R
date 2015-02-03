
library("basefun")
n <- 100
set.seed(29)

### one-dim problem
## set-up simple regression problem
x <- sort(runif(100))
y <- x^2 + rnorm(length(x), sd = .1) 
## set-up monotone Bernstein polynom basis
Bb <- Bernstein_basis(order = 5, ui = "increasing", varname = "x")
## evaluate basis
X1 <- model.matrix(Bb, data = data.frame(x = x))
X2 <- Bb(x)
stopifnot(all.equal(X1, X2))
## fit model
m1 <- lm(y ~  X1 - 1)
m2 <- lm(y ~  Bb(x) - 1, data = data.frame(y = y, x = x))
Bb0 <- Bernstein_basis(order = 5, ui = "zeroint", varname = "x")
X0 <- model.matrix(Bb0, data = data.frame(x = x))
m0 <- lm(y ~  X0)
stopifnot(all.equal(fitted(m1), fitted(m2)))
stopifnot(all.equal(fitted(m1), fitted(m0)))

stopifnot(all.equal(coef(m1), coef(m2), check.attributes = FALSE))
## generate new data from support of basis
xn <- generate(Bb, n = 100)
## compute estimated regression function
p1 <- predict(Bb, newdata = data.frame(x = xn), coef = coef(m1))
p2 <- predict(m2, newdata = data.frame(x = xn)) 
stopifnot(all.equal(c(p1), p2, check.attributes = FALSE))
## compute derivative of estimated regression function
dp1 <- predict(Bb, newdata = data.frame(x = xn), coef = coef(m1), deriv = 1)
dp0 <- predict(Bb0, newdata = data.frame(x = xn), coef = coef(m0)[-1], deriv = 1)
stopifnot(all.equal(dp1, dp0))
dp12 <- predict(Bb, newdata = data.frame(x = xn), coef = coef(m1), deriv = 1, integrate = TRUE)
unique(c(p1 - dp12)) ### the same up to a constant

### two-dim (ANCOVA) problem
gf <- gl(3, 1)
g <- sample(gf, length(x), replace = TRUE)
y <- x^2 + (g == "2") * sin(x) + (g == "3") * cos(x) + rnorm(length(x), sd = .05)
## generate a basis for a factor (treatment contrasts)
## this is equal to model.matrix(~ gf)
gb <- as.basis(~ g, remove_intercept = FALSE, data = data.frame(g = gf))
## join the two bases by the kronecker product
bb <- b(b1 = Bb, b2 = gb)
## evaluate new two-dim basis
X1 <- model.matrix(bb, data = data.frame(x = x, g = g))
X2 <- bb(data.frame(x = x, g = g))
stopifnot(all.equal(X1, X2))
## fit model
m1 <- lm(y ~  X1 - 1)
m2 <- lm(y ~  bb(data.frame(x = x, g = g)) - 1, data = data.frame(y = y, x = x, g = g))
stopifnot(all.equal(coef(m1), coef(m2), check.attributes = FALSE))
## compute estimated regression functions
d <- generate(bb, n = 100)
## for each group
p1 <- sapply(gf, function(l) predict(bb, newdata = data.frame(x = d$x, g = l), coef = coef(m1)))
## the same via _linear array_ approach
p2 <- predict(bb, newdata = d, coef(m1))
## brute force; 2 times
p3 <- matrix(predict(bb, newdata = do.call(expand.grid, d), coef(m1)), ncol = nlevels(gf))
p4 <- matrix(predict(m2, newdata = do.call(expand.grid, d)), ncol = nlevels(gf))
stopifnot(all.equal(p1, p2))
stopifnot(all.equal(p2, p3))
stopifnot(all.equal(p3, p4))
## compute derivative wrt the first element
dp2 <- predict(bb, newdata = d, coef(m1), deriv = c(x = 1))
