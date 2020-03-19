# Tests for COLR models

## Depencies
## IGNORE_RDIFF_BEGIN
library("tramnet")
library("survival")
library("sandwich")
## IGNORE_RDIFF_END
options(digits = 3)

## Exact and Right censored
data("GBSG2", package = "TH.data")
GBSG2$surv <- with(GBSG2, Surv(time, cens))
x <- matrix(1 * (GBSG2$horTh == "yes"), ncol = 1)
colnames(x) <- "horTh"

yCOLR <- Colr(surv ~ 1, data = GBSG2, log_first = TRUE, order = 10)
modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
yCOLRb <- Colr(surv ~ horTh, data = GBSG2, log_first = TRUE, order = 10)
max(abs(coef(yCOLRb, with_baseline = FALSE) -
        coef(modCOLR, with_baseline = FALSE)))
logLik(yCOLRb)
logLik(modCOLR)
-modCOLR$result$value
logLik(modCOLR, newdata = tramnet:::.get_tramnet_data(modCOLR)[1:10, ])

## methods
coef(modCOLR, tol = 0, with_baseline = TRUE)
logLik(modCOLR)
resid(modCOLR)[1:10]
predict(modCOLR, type = "distribution", q = 1)[, 1:10]
predict(modCOLR, type = "quantile", prob = 0.5)[, 1:10]
head(simulate(modCOLR))
head(estfun(modCOLR))
plot(modCOLR, type = "survivor")
plot(modCOLR, type = "density", K = 120)
print(modCOLR)


if (FALSE) {
  ## left censored
  GBSG2$cens <- as.integer(GBSG2$cens)
  GBSG2$cens[GBSG2$time < 200] <- 2
  GBSG2$time[GBSG2$cens == 2] <- 100

  yCOLR <- Colr(Surv(time, time, cens, type = "interval") ~ 1, data = GBSG2, log_first = TRUE, order = 10)
  modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
  yCOLRb <- Colr(Surv(time, time, cens, type = "interval") ~ horTh, data = GBSG2, log_first = TRUE, order = 10)
  max(abs(coef(yCOLRb, with_baseline = FALSE) -
          coef(modCOLR, with_baseline = FALSE)))
  logLik(yCOLRb)
  logLik(modCOLR)

  ## methods
  coef(modCOLR, tol = 0, with_baseline = TRUE)
  logLik(modCOLR)
  resid(modCOLR)[1:10]
  predict(modCOLR, type = "distribution", q = 1)[, 1:10]
  predict(modCOLR, type = "quantile", prob = 0.5)[, 1:10]
  head(simulate(modCOLR))
  head(estfun(modCOLR))
  plot(modCOLR, type = "survivor")
  plot(modCOLR, type = "density", K = 120)
  print(modCOLR)

  ## Unconditional, stratified
  yCOLR <- Colr(surv | horTh ~ 1, data = GBSG2)
  modCOLR <- tramnet(yCOLR, x = matrix(0, nrow = nrow(GBSG2)), lambda = 0, alpha = 0)
  logLik(yCOLR)
  logLik(modCOLR)
  resid(modCOLR)[1:10]
  predict(modCOLR, type = "distribution", q = 1)[, 1:10]
  predict(modCOLR, type = "quantile", prob = 0.5)[, 1:10]
  head(simulate(modCOLR))
  head(estfun(modCOLR))
  plot(modCOLR, type = "survivor")
  plot(modCOLR, type = "density", K = 120)
  print(modCOLR)


  ## interval censored
  GBSG2$time2 <- GBSG2$time + 50
  GBSG2$cens[which(GBSG2$cens == 1)[1:100]] <- 3

  yCOLR <- Colr(Surv(time, time2, cens, type = "interval") ~ 1, data = GBSG2, log_first = TRUE, order = 10)
  modCOLR <- tramnet(yCOLR, x, lambda = 0, alpha = 0)
  yCOLRb <- Colr(Surv(time, time2, cens, type = "interval") ~ horTh, data = GBSG2, log_first = TRUE, order = 10)
  max(abs(coef(yCOLRb, with_baseline = FALSE) -
            coef(modCOLR, with_baseline = FALSE)))
  logLik(yCOLRb)
  logLik(modCOLR)
}
