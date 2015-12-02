	
library("sltm")
library("multcomp")

(a <- sltm(Sepal.Length | Species ~ Sepal.Width, data = iris, trafo = "log"))
summary(a)
coef(a)
vcov(a)
predict(a)

library("survival")
data("GBSG2", package = "TH.data")
storage.mode(GBSG2$time) <- "double"

(b <- sltm(time ~ horTh , data = GBSG2))
summary(b)
coef(b)
vcov(b)
predict(b)
predict(b, newdata = data.frame(horTh = unique(GBSG2$horTh)), q = 500:510, type = "trafo")

(b <- sltm(time ~ 1, data = GBSG2))
summary(b)
coef(b)
vcov(b)
predict(b)
predict(b, q = 500:510, type = "distribution")

(b <- sltm(time | tgrade ~ 1, data = GBSG2))
summary(b)
coef(b)
vcov(b)
predict(b)
predict(b, newdata = data.frame(tgrade = unique(GBSG2$tgrade)), q = 500:510, type = "distribution")

(cc <- sltm(Surv(time, cens) | tgrade ~ horTh + menostat + pnodes, data = GBSG2, method = "cloglog"))
summary(cc)
coef(cc)
vcov(cc)
predict(cc)
predict(cc, newdata = expand.grid(horTh = sort(unique(GBSG2$horTh)), 
                                  menostat = sort(unique(GBSG2$menostat))[1L],
                                  pnodes = 30L, 
                                  tgrade = sort(unique(GBSG2$tgrade))[2L]), 
         q = 500:501, type = "trafo")

(d <- sltm(Surv(time, cens) ~ 1, data = GBSG2, method = "cloglog"))
summary(d)
coef(d)
vcov(d)
predict(d)

cb <- confband(d, newdata = data.frame(1), calpha = univariate_calpha(), type = "surv")

layout(matrix(1:2, ncol = 2))
plot(survfit(Surv(time, cens) ~ 1, data = GBSG2))
plot(d, newdata = data.frame(1), type = "surv")
lines(cb[,"q"], cb[,"lwr"])
lines(cb[,"q"], cb[,"upr"])
