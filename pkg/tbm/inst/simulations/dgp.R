
library("tram")

set.seed(27031105)

n <- 1000
order <- 6
sup <- c(-4, 6)
bds <- c(-8, 8)
N <- 500
beta <- c(-2, 2, -1)

x1 <- gl(2, 1)
x2 <- seq(from = 0, to = 1, length.out = n)
d <- expand.grid(x1 = x1, x2 = x2)
d$int <- 1
d$f2 <- with(d, sin(x2 * 2 * pi) * (1 + x2))
X <- model.matrix(~ x1 * x2, data = d)[,-1]
d$y <- rnorm(nrow(d), mean = X %*% beta)

m0 <- BoxCox(y ~ 1, data = d, order = order, 
             model_only = TRUE, support = sup,
             bounds = bds)
cf0 <- seq(from = sup[1], to = sup[2], length = order + 1) + 
       sin(seq(from = 0, to = 2*pi, length = order + 1)) * (seq(from = 0, to = 2*pi, length = order + 1) <= pi) * 2

m1 <- BoxCox(y ~ x1 * x2, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m1)
cf[] <- c(cf0, beta)
coef(m1) <- cf

m2 <- BoxCox(y ~ x1 * f2, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m2)
cf[] <- c(cf0, beta)
coef(m2) <- cf

m3 <- BoxCox(y | x1 * x2 ~ 1, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m3)

cf1 <- sin(seq(from = 0, to = pi / 2, length.out = order + 1)) * beta[1]
cf2 <- sqrt(seq(from = 0, to = 2, length.out = order + 1)) / sqrt(2) * beta[2]
cf21 <- sin(seq(from = 0, to = pi / 2, length.out = order + 1)) * beta[3]

cf[] <- c(cf0, cf1, cf2, cf21)
coef(m3) <- cf

m4 <- BoxCox(y | x1 * f2 ~ 1, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m4)

cf[] <- c(cf0, cf1, cf2, cf21)
coef(m4) <- cf

q <- mkgrid(m1, n = N)[["y"]]
x2 <- 0:10 / 10
nd <- expand.grid(x1 = x1, x2 = x2)
ndq <- expand.grid(y = q, x1 = x1, x2 = x2)
ndq$xconfig <- with(ndq, interaction(x1, factor(x2)))
nd$f2 <- with(nd, sin(x2 * 2 * pi) * (1 + x2))

ndq$d1 <- c(predict(m1, newdata = nd, type = "density", q = q))
ndq$d2 <- c(predict(m2, newdata = nd, type = "density", q = q))
ndq$d3 <- c(predict(m3, newdata = nd, type = "density", q = q))
ndq$d4 <- c(predict(m4, newdata = nd, type = "density", q = q))

ndq$t1 <- c(predict(m1, newdata = nd, type = "trafo", q = q))
ndq$t2 <- c(predict(m2, newdata = nd, type = "trafo", q = q))
ndq$t3 <- c(predict(m3, newdata = nd, type = "trafo", q = q))
ndq$t4 <- c(predict(m4, newdata = nd, type = "trafo", q = q))

col <- rep(grey.colors(length(x2), alpha = .9), each = 2)
xlim <- bds

save(xlim, col, ndq, file = "plot_dgp.rda")

library("lattice")
pdf("dgp.pdf")
xyplot(d1 ~ y | x1, data = ndq, type = "l", groups = xconfig, col = col, lwd =
2, xlim = xlim)
xyplot(d2 ~ y | x1, data = ndq, type = "l", groups = xconfig, col = col, lwd =
2, xlim = xlim)
xyplot(d3 ~ y | x1, data = ndq, type = "l", groups = xconfig, col = col, lwd =
2, xlim = xlim)
xyplot(d4 ~ y | x1, data = ndq, type = "l", groups = xconfig, col = col, lwd =
2, xlim = xlim)
dev.off()

y1 <- simulate(m1, newdata = d, n = 300)
y2 <- simulate(m2, newdata = d, n = 300)
y3 <- simulate(m3, newdata = d, n = 300)
y4 <- simulate(m4, newdata = d, n = 300)

y <- list(y1 = y1, y2 = y2, y3 = y3, y4 = y4)
model <- list(m1 = m1, m2 = m2, m3 = m3, m4 = m4)

save(order, sup, bds, d, y, cf0, cf1, cf2, cf21, model, file = "dgp.rda")

