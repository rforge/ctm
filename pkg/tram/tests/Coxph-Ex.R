
library("tram")
library("survival")
library("trtf")

data("GBSG2", package = "TH.data")

cmod <- coxph(Surv(time, cens) ~ progrec + pnodes + strata(horTh, tgrade),
               data = GBSG2)
summary(cmod)

Cmod <- Coxph(Surv(time, cens) | 0 + horTh:tgrade ~ progrec + pnodes, 
              data = GBSG2)
summary(Cmod)

coef(cmod)
coef(Cmod)

vcov(cmod)
vcov(Cmod)

cmod <- coxph(Surv(time, cens) ~ ., data = GBSG2)
Cmod <- Coxph(Surv(time, cens) ~ ., data = GBSG2)

coef(cmod)
coef(Cmod)

diag(vcov(cmod))
diag(vcov(Cmod))

if (FALSE) {
fm <- Coxph(Surv(time, cens) ~ 1, data = GBSG2, asFamily = TRUE)
bmod <- glmboost(Surv(time, cens) ~ ., data = GBSG2, family = fm, 
                 control = boost_control(nu = .75, mstop = 250, trace = TRUE), 
                 center = FALSE)
cf <- coef(bmod)
ns <- nuisance(bmod)

cf2 <- c(ns[names(ns) != "(Intercept)"], -cf)

logLik(as.mlt(Cmod), parm = cf2)
logLik(Cmod)
coef(as.mlt(Cmod))
cf2

risk(bmod)

fm <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2, asFamily = TRUE)
bmod <- glmboost(Surv(time, cens) ~ . - horTh, data = GBSG2, family = fm, 
                 control = boost_control(nu = .75, mstop = 250, trace = TRUE), center = FALSE)
cf <- coef(bmod)
ns <- nuisance(bmod)

cf2 <- c(ns[names(ns) != "(Intercept)"], -cf)

logLik(as.mlt(Cmod), parm = cf2)
logLik(Cmod)
coef(as.mlt(Cmod))
cf2

risk(bmod)

}

cmod <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2)
(tmod <- trafotree(cmod, formula = Surv(time, cens) ~ horTh | ., data = GBSG2))
logLik(tmod)
