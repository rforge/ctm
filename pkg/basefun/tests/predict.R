
library("basefun")

x <- 1:5
y <- as.double(1:4)
g <- gl(3, 1)
d <- expand.grid(x = x, y = y, g = g)

cb <- c(logx = log_basis(x, varname = "x"), 
        X = as.basis(~ y + g, data = expand.grid(y = y, g = g)))

X <- model.matrix(cb, data = d)
stopifnot(nrow(X) == nrow(d))

p <- predict(cb, newdata = d, coef = rep(1, ncol(X)))
stopifnot(length(p) == nrow(d))

(p2 <- predict(cb, newdata = list(x = x, y = y, g = g), coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p2), check.attributes = FALSE))

(p3 <- predict(cb, newdata = list(logx = list(x = x), X = expand.grid(y = y, g = g)), 
               coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p3), check.attributes = FALSE))

(p4 <- predict(cb, newdata = generate(cb, 4), coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p4), check.attributes = FALSE))

XX <- model.matrix(cb[["X"]], data = list(y = y, g = g))

pX <- predict(cb[["X"]], newdata = list(y = y, g = g), coef = rep(1, ncol(X) - 1))

pX2 <- predict(cb[["X"]], newdata = expand.grid(y = y, g = g), coef = rep(1, ncol(X) - 1))

stopifnot(all.equal(pX, c(pX2), check.attributes = FALSE))

bb <- b(logx = log_basis(x, varname = "x"),
        X = as.basis(~ y + g, data = expand.grid(y = y, g = g)))

X <- model.matrix(bb, data = d)
stopifnot(nrow(X) == nrow(d))  

p <- predict(bb, newdata = d, coef = rep(1, ncol(X)))
stopifnot(length(p) == nrow(d))

(p2 <- predict(bb, newdata = list(x = x, y = y, g = g), coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p2), check.attributes = FALSE))
