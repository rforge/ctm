
library("ltm")
library("multcomp")
library("survival")
data("GBSG2", package = "TH.data")
storage.mode(GBSG2$time) <- "double"

mltm <- ltm(Surv(time, cens) | tgrade ~ horTh + menostat + pnodes, data = GBSG2, method = "cloglog",
            negative_lp = FALSE, order = 5)
mcox <- coxph(Surv(time, cens) ~ horTh + menostat + pnodes + strata(tgrade), data = GBSG2)

print(mltm)
print(mcox)

summary(mltm)
summary(mcox)

confint(mltm, calpha = univariate_calpha())
confint(mcox)

coef(mltm)
coef(mcox)

vcov(mltm)
vcov(mcox)

summary(predict(mltm, newdata = GBSG2, type = "lp") - model.matrix(mcox) %*% coef(mcox))

nd <- expand.grid(horTh = sort(unique(GBSG2$horTh)), 
                  menostat = sort(unique(GBSG2$menostat))[1L],
                  pnodes = 30L, 
                  tgrade = sort(unique(GBSG2$tgrade))[2L])

plot(survfit(mcox, newdata = nd, col = "grey"))
plot(mltm, newdata = nd, type = "surv", add = TRUE, col = "red")

