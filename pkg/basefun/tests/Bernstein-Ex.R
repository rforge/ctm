
library("basefun")

### check the approximation of a number of functions
f1 <- function(x) qnorm(pchisq(x, df = 3))
fun <- list(sin, cos, sqrt, log, f1)
dfun <- list(cos, sin, function(x) -.5 * x^(-3/2), function(x) 1/x, 
             function(x) dnorm(pchisq(x, df = 3)) * dchisq(x, df = 3))
ord <- 3:10
x <- seq(from = 0.01, to = 2*pi - 0.01, length.out = 100)
for (i in 1:length(fun)) {
    for (o in ord) {
        y <- fun[[i]](x)
        Bb <- Bernstein_basis(order = o, varname = "x", support = c(0, 2*pi))
        m <- lm(y ~ Bb(x) - 1, data = data.frame(y = y, x = x))
        R <- summary(m)$r.squared
        layout(matrix(1:2, ncol = 2))
        plot(x, fun[[i]](x), type = "l", col = "red", main = paste(deparse(fun[[i]]), o, R, sep = ":"))
        lines(x, fitted(m))
        plot(x, dfun[[i]](x), type = "l", col = "red", main = paste(deparse(fun[[i]]), o, R, sep = ":"))
        lines(x, predict(Bb, newdata = data.frame(x = x), deriv = 1, coef = coef(m)))
    }
}

