
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
o <- 6
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

o <- 6
Bs <- Bernstein_basis(order = o, var = "waiting", ui = "incre",
                      support = range(faithful$waiting) + c(-5, 5))
m <- model(Bs)
mod <- mlt(m, data = faithful)

library("multcomp")

mp <- parm(coef(mod), vcov(mod))
y <- generate(mod, 50)$waiting
mc <- confint(glht(mp, linfct = model.matrix(mod$model, 
    data = data.frame(waiting = y))))
umc <- confint(glht(mp, linfct = model.matrix(mod$model, 
    data = data.frame(waiting = y))), calpha = qnorm(.975))
p <- mod$model$todistr$p
plot(y, p(mc$confint[, "Estimate"]), type = "l")
lines(y, p(mc$confint[, "lwr"]))
lines(y, p(mc$confint[, "upr"]))
lines(y, p(umc$confint[, "lwr"]))
lines(y, p(umc$confint[, "upr"]))

library("survival")
cm <- coxph(Surv(waiting, rep(TRUE, nrow(faithful))) ~ 1, data = faithful)
plot(survfit(cm))
lines(y, 1 - p(mc$confint[, "Estimate"]), col = "red")
lines(y, 1 - p(mc$confint[, "lwr"]), col = "red")
lines(y, 1 - p(mc$confint[, "upr"]), col = "red")
