
library("ltm")
library("multcomp")

(a <- ltm(Sepal.Length | Species ~ Sepal.Width, data = iris, trafo = "log"))
summary(a)
coef(a)
vcov(a)
predict(a)

library("survival")
data("GBSG2", package = "TH.data")
storage.mode(GBSG2$time) <- "double"

(b <- ltm(time ~ horTh , data = GBSG2))
summary(b)
coef(b)
vcov(b)
predict(b)
predict(b, newdata = data.frame(horTh = unique(GBSG2$horTh)), q = 500:510, type = "trafo")

(b <- ltm(time ~ 1, data = GBSG2))
summary(b)
coef(b)
vcov(b)
predict(b)
predict(b, q = 500:510, type = "distribution")

(b <- ltm(time | tgrade ~ 1, data = GBSG2))
summary(b)
coef(b)
vcov(b)
predict(b)
predict(b, newdata = data.frame(tgrade = unique(GBSG2$tgrade)), q = 500:510, type = "distribution")

(cc <- ltm(Surv(time, cens) | tgrade ~ horTh + menostat + pnodes, data = GBSG2, method = "cloglog"))
summary(cc)
coef(cc)
vcov(cc)
predict(cc)
predict(cc, newdata = expand.grid(horTh = sort(unique(GBSG2$horTh)), 
                                  menostat = sort(unique(GBSG2$menostat))[1L],
                                  pnodes = 30L, 
                                  tgrade = sort(unique(GBSG2$tgrade))[2L]), 
         q = 500:501, type = "trafo")

(d <- ltm(Surv(time, cens) ~ 1, data = GBSG2, method = "cloglog"))
summary(d)
coef(d)
vcov(d)
predict(d)

confint(d, calpha = univariate_calpha())

d <- as.mlt(d)
cf <- coef(d)
v <- vcov(d)

prm <- parm(cf, v)
K <- diag(length(cf))
rownames(K) <- colnames(K) <- names(cf)
ci <- confint(glht(prm, linfct = K), calpha = qnorm(.975))

lwr <- d
upr <- d
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

