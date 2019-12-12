# Tests for Lehmann alternative models

## Dependencies
library("survival")
library("tramnet")
options(digits = 4)

## Data
data("GBSG2", package = "TH.data")
GBSG2$surv <- with(GBSG2, Surv(time, rep(1, nrow(GBSG2))))
x <- matrix(1 * (GBSG2$horTh == "yes"), ncol = 1) ## create matrix of covariates
colnames(x) <- "horTh"

## Exact
yLehmann <- Lehmann(surv ~ 1, data = GBSG2, log_first = TRUE, order = 10)
modLehmann <- tramnet(yLehmann, x, lambda = 0, alpha = 0)

yLehmannb <- Lehmann(surv ~ horTh, data = GBSG2, log_first = TRUE, order = 10)
max(abs(coef(yLehmannb, with_baseline = FALSE) -
        coef(modLehmann, with_baseline = FALSE)))
logLik(yLehmannb)
logLik(modLehmann)

## left censored
GBSG2$cens <- 1
GBSG2$cens <- as.integer(GBSG2$cens)
GBSG2$cens[GBSG2$time < 100] <- 2
GBSG2$time[GBSG2$cens == 2] <- 100

yCOLR <- Lehmann(Surv(time, time, cens, type = "interval") ~ 1, data = GBSG2, log_first = TRUE, order = 10)
modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
yCOLRb <- Lehmann(Surv(time, time, cens, type = "interval") ~ horTh, data = GBSG2, log_first = TRUE, order = 10)
max(abs(coef(yCOLRb, with_baseline = FALSE) -
        coef(modCOLR, with_baseline = FALSE)))
logLik(yCOLRb)
logLik(modCOLR)

