
library("tram")
library("survival")

data("GBSG2", package = "TH.data")

fm <- Surv(time, cens) ~ pnodes + age
tfm <- Surv(time, cens) ~ pnodes + age
fms <- Surv(time, cens) ~ pnodes + age + strata(horTh)
tfms <- Surv(time, cens) | 0 + horTh ~ pnodes + age

(smod <- survreg(fm, data = GBSG2, dist = "weibull"))
(Smod <- Survreg(tfm, data = GBSG2, dist = "weibull"))
Smod$invscale

summary(Smod)

(smod <- survreg(fms, data = GBSG2, dist = "weibull"))
(Smod <- Survreg(tfms, data = GBSG2, dist = "weibull"))
Smod$invscale

(smod <- survreg(fm, data = GBSG2, dist = "exponential"))
(Smod <- Survreg(tfm, data = GBSG2, dist = "exponential"))
Smod$invscale

try(smod <- survreg(fms, data = GBSG2, dist = "exponential"))
(Smod <- Survreg(tfms, data = GBSG2, dist = "exponential"))
Smod$invscale

(smod <- survreg(fm, data = GBSG2, dist = "rayleigh"))
(Smod <- Survreg(tfm, data = GBSG2, dist = "rayleigh"))
Smod$invscale

try(smod <- survreg(fms, data = GBSG2, dist = "rayleigh"))
(Smod <- Survreg(tfms, data = GBSG2, dist = "rayleigh"))
Smod$invscale

(smod <- survreg(fm, data = GBSG2, dist = "lognormal"))
(Smod <- Survreg(tfm, data = GBSG2, dist = "lognormal"))
Smod$invscale

(smod <- survreg(fms, data = GBSG2, dist = "lognormal"))
(Smod <- Survreg(tfms, data = GBSG2, dist = "lognormal"))
Smod$invscale

(smod <- survreg(fm, data = GBSG2, dist = "loglogistic"))
(Smod <- Survreg(tfm, data = GBSG2, dist = "loglogistic"))
Smod$invscale

(smod <- survreg(fms, data = GBSG2, dist = "loglogistic"))
(Smod <- Survreg(tfms, data = GBSG2, dist = "loglogistic"))
coef(Smod) * Smod$invscale
Smod$invscale

(tobinfit <- survreg(Surv(durable, durable>0, type='left') ~ age + quant,
                         data=tobin, dist='gaussian'))


(tobinfit2 <- Survreg(Surv(durable, durable>0, type='left') ~ age + quant,
                         data=tobin, dist='gaussian'))

coef(tobinfit2) * tobinfit2$invscale
