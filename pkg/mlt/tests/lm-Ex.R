
library("mlt")

n <- 1000	
y <- rnorm(n, 2, 1.5)
d <- data.frame(y = y)

cnst <- list(ui = Diagonal(2), ci = c(-Inf, 0))
m <- model(as.bases(~ 1 + y, vars = d, constraint = cnst))
m(d)
m(d, bresponse = list(deriv = 1))

o <- mlt(m, data = d)

s2ml <- sqrt(var(y) * (n - 1) / n)
1 / coef(o)[2] - s2ml

-coef(o)[1] / coef(o)[2] - mean(y)

x <- runif(n, max = 2 * pi)
y <- rnorm(n, sin(x), .25)
d <- data.frame(y = y, x = x)

plot(x, y)

m <- model(as.bases(~ y - 1, vars = d, constraint = list(ui = Diagonal(1), ci = 0)), 
           shift = Bernstein_basis(order = 10, support = c(0, 2*pi), var = "x"))

o <- mlt(m, data = d)
1 / coef(o)[1]
p <- predict(Bernstein_basis(order = 10, support = c(0, 2*pi), var = "x"), 
             newdata = data.frame(x = x), coef = coef(o)[-1])
plot(x, y)
lines(sort(x), -p[order(x)] / coef(o)[1], lwd = 2, col = "red")


x <- runif(n, max = 2 * pi)
y <- rnorm(n, 2, 1.1 + sin(x) / 2)
d <- data.frame(y = y, x = x)

plot(x, y)

m <- model(as.bases(~ 1 + y, vars = d, constraint = cnst), 
           Bernstein_basis(order = 10, support = c(0, 2*pi), var = "x"))

o <- mlt(m, data = d)

