# Tests for Coxph models

## Dependencies
## IGNORE_RDIFF_BEGIN
library("tramnet")
library("survival")
## IGNORE_RDIFF_END
options(digits = 3)

## Exact and Right censored
data("GBSG2", package = "TH.data")
GBSG2$surv <- with(GBSG2, Surv(time, cens))
x <- matrix(1 * (GBSG2$horTh == "yes"), ncol = 1) ## create matrix of covariates
colnames(x) <- "horTh"

yCOLR <- Coxph(Surv(time, cens) ~ 1, data = GBSG2, log_first = TRUE, order = 10)
modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
yCOLRb <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2, log_first = TRUE, order = 10)
max(abs(coef(yCOLRb, with_baseline = FALSE) -
        coef(modCOLR, with_baseline = FALSE)))
logLik(yCOLRb)
logLik(modCOLR)


if (FALSE) {
  ## left censored
  GBSG2$cens <- as.integer(GBSG2$cens)
  GBSG2$cens[GBSG2$time < 100] <- 2
  GBSG2$time[GBSG2$cens == 2] <- 100

  yCOLR <- Coxph(Surv(time, time, cens, type = "interval") ~ 1, data = GBSG2,
                 log_first = TRUE, order = 10)
  modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
  yCOLRb <- Coxph(Surv(time, time, cens, type = "interval") ~ horTh, data = GBSG2,
                  log_first = TRUE, order = 10)
  max(abs(coef(yCOLRb, with_baseline = FALSE) -
            coef(modCOLR, with_baseline = FALSE)))
  logLik(yCOLRb)
  logLik(modCOLR)
}

if (FALSE) {
  ## interval censored
  GBSG2$time2 <- GBSG2$time + 50
  GBSG2$cens[which(GBSG2$cens == 1)[1:100]] <- 3

  yCOLR <- Coxph(Surv(time, time2, cens, type = "interval") ~ 1, data = GBSG2,
                 log_first = TRUE, order = 10)
  modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
  yCOLRb <- Coxph(Surv(time, time2, cens, type = "interval") ~ horTh, data = GBSG2,
                  log_first = TRUE, order = 10)
  max(abs(coef(yCOLRb, with_baseline = FALSE) -
            coef(modCOLR, with_baseline = FALSE)))
  logLik(yCOLRb)
  logLik(modCOLR)
}
