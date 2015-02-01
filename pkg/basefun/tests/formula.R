
library("basefun")

n <- 100
x <- runif(n)
y <- runif(n)
g <- gl(4, n / 4)
d <- data.frame(y = y, x = x, g = g)
de <- d[-(1:nrow(d)),]

b1 <- as.bases(~ x + g, data = de)
b2 <- Bernstein_basis(order = 4, var = "y", ui = "incre")

b12 <- b(b1 = b1, b2 = b2)
c12 <- c(b1 = b1, b2 = b2)

dim(b12(d))
dim(c12(d))

tmp <- c(b12 = b12, c12 = c12)
dim(tmp(d))
