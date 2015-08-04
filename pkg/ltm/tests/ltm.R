
library("ltm")
library("multcomp")

coef(a <- ltm(Sepal.Length ~ Sepal.Width | Species, data = iris, trafo = "log"))
predict(a)

library("survival")
data("GBSG2", package = "TH.data")
storage.mode(GBSG2$time) <- "double"

b <- ltm(time ~ horTh , data = GBSG2)

predict(b)

predict(b, newdata = data.frame(horTh = unique(GBSG2$horTh)), q = 100:110, type = "trafo")

b <- ltm(time ~ 1, data = GBSG2)

predict(b, q = 100:110, type = "distribution")

b <- ltm(time ~ 1 | tgrade, data = GBSG2)

predict(b, newdata = data.frame(tgrade = unique(GBSG2$tgrade)), q = 100:110, type = "distribution")

cc <- ltm(Surv(time, cens) ~ horTh + menostat + pnodes | tgrade, data = GBSG2, method = "cloglog")

predict(cc, newdata = expand.grid(horTh = c("no", "yes"), menostat = "Pre",
                           pnodes = 100, tgrade = "II"), q = 100:101, type = "trafo")

d <- ltm(Surv(time, cens) ~ 1, data = GBSG2, method = "cloglog")

confint(d, calpha = univariate_calpha())

class(d) <- class(d)[-1]
cf <- coef(d)
v <- vcov(d)


prm <- parm(cf, v)
K <- diag(length(cf))
rownames(K) <- colnames(K) <- names(cf)
ci <- confint(glht(prm, linfct = K), calpha = qnorm(.975))

lwr <- d
class(lwr) <- class(lwr)[-1]
upr <- d
class(upr) <- class(upr)[-1]
coef(lwr) <- ci$confint[, "lwr"]
coef(upr) <- ci$confint[, "upr"]

tm <- 10:2700

s <- seq(from = min(GBSG2$time), to = max(GBSG2$time), length = length(cf))
plot(tm, predict(d, q = tm), type = "l", ylim = range(ci$confint))
points(s, cf, col = "red")

lines(tm, predict(lwr, q = tm))
points(s, coef(lwr), col = "blue")

lines(tm, predict(upr, q = tm))
points(s, coef(upr), col = "green")

plot(tm, 1 - predict(d, q = tm, type = "distr"), type = "l")

lines(tm, 1 - predict(lwr, q = tm, type = "distr"))

lines(tm, 1 - predict(upr, q = tm, type = "distr"))



plot(survfit(Surv(time, cens) ~ 1, data = GBSG2))
lines(tm, 1 - predict(d, q = tm, type = "distr"))
lines(tm, 1 - predict(lwr, q = tm, type = "distr"))
lines(tm, 1 - predict(upr, q = tm, type = "distr"))

