
library("tram")
library("survival")
library("trtf")

### Windows diffs...
options(digits = 3)

data("GBSG2", package = "TH.data")

cmod <- coxph(Surv(time, cens) ~ progrec + pnodes + strata(horTh, tgrade),
               data = GBSG2)
summary(cmod)

Cmod <- Coxph(Surv(time, cens) | 0 + horTh:tgrade ~ progrec + pnodes, 
              data = GBSG2)
summary(Cmod)

Cmod_lf <- Coxph(Surv(time, cens) | 0 + horTh:tgrade ~ progrec + pnodes, 
                data = GBSG2, log_first = TRUE)
summary(Cmod_lf)

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

cmod <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2)
(tmod <- trafotree(cmod, formula = Surv(time, cens) ~ horTh | ., data = GBSG2))
logLik(tmod)
