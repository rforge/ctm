
library("tram")
library("survival")
library("mboost")
library("trtf")

data("GBSG2", package = "TH.data")

(cmod <- coxph(Surv(time, cens) ~ age + pnodes + strata(horTh, tgrade),
               data = GBSG2))

summary(cmod)

(Cmod <- Coxph(Surv(time, cens) | 0 + horTh:tgrade ~ age + pnodes, 
              data = GBSG2))
summary(Cmod)

coef(Cmod)
vcov(Cmod)
nobs(Cmod)

(cmod <- coxph(Surv(time, cens) ~ ., data = GBSG2))

fm <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2, asFamily = TRUE)
bmod <- glmboost(Surv(time, cens) ~ . - horTh, data = GBSG2, family = fm, 
                 control = boost_control(nu = .5), center = FALSE)
coef(bmod)
nuisance(bmod)
risk(bmod)

cmod <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2)
(tmod <- trafotree(cmod, formula = Surv(time, cens) ~ horTh | ., data = GBSG2))

logLik(tmod)
