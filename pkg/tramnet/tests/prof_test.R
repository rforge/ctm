# Test for profiling cfx trajectory functions

## IGNORE_RDIFF_BEGIN
library("tramnet")
library("penalized")
library("survival")
## IGNORE_RDIFF_END
options(digits = 3)

## --- Comparison with penalized
data("nki70", package = "penalized")
nki70$resp <- with(nki70, Surv(time, event))
x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
fit2 <- tramnet(y, x, lambda = 0, alpha = 1)
(pfl <- prof_lambda(fit2))
plot_path(pfl)
fit3 <- tramnet(y, x, lambda = 1, alpha = 1)
(pfa <- prof_alpha(fit3))
plot_path(pfa)
