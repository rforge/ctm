
library("tramnet")
library("penalized")
library("survival")
options(digits = 3)

## --- Comparison with penalized
data("nki70", package = "penalized")
resp <- with(nki70, Surv(time, event))
x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
fit <- penalized(response = resp, penalized = x, lambda1 = 1, lambda2 = 0,
                 standardize = FALSE, data = nki70)
y <- Coxph(Surv(time, event) ~ 1, data = nki70, order = 10, log_first = TRUE)
fit2 <- tramnet(y, x, lambda = 1, alpha = 1) ## L1 only
coef(fit)
coef(fit2)


