
library("mlt")
library("lattice")

data("faithful")

aic <- numeric(20)

for (o in 2:(length(aic) + 1)) {
Bs <- Bernstein_basis(order = o, var = "waiting", ui = "incre", 
                      support = range(faithful$waiting) + c(-5, 5))
m <- model(Bs)
mod <- mlt(m, data = faithful)

aic[o - 1] <- AIC(mod)

plot(mod, ~ waiting, what = "prob",  
     col = "red", pch = 21, main = paste("order", o, "aic", aic[o - 1]))
lines(ecdf(faithful$waiting))

}

plot(aic)

o <- which.min(aic) + 1
Bs <- Bernstein_basis(order = o, var = "waiting", ui = "incre",
                      support = range(faithful$waiting) + c(-5, 5))
m <- model(Bs)
mod <- mlt(m, data = faithful)

abline(h = AIC(mod))

p <- predict(mod, newdata = data.frame(1))

y <- generate(mod, 100)$waiting

tr <- p[[1]](y)
tr1 <- p[[1]](y, deriv = c(waiting = 1))
plot(y, mod$todistr$d(tr) * tr1, type = "l", col = "red", lwd = 3)
lines(density(faithful$waiting))
lines(rug(faithful$waiting))
abline(h = 0)

d <- predict(mod, newdata = data.frame(1))

y <- generate(mod, 100)$waiting

tr <- d[[1]](y, type = "density")
lines(y, tr, type = "l", col = "blue", lwd = 1)

p <- 1:99 / 100

plot(p, d[[1]](p, type = "quantile"))
lines(p, quantile(faithful$waiting, p))
