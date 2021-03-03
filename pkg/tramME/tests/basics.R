## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")

## -- Model setup with initpar
ip <- list(beta = c(1.5, 0.08, 0.2))
mod <- LmME(dist ~ speed, data = cars, initpar = ip, nofit = TRUE)
## NOTE: initpars are not set as actual model parameters...
chkerr(chkeq(ip$beta, coef(mod, with_baseline = TRUE), check.attributes = FALSE))
## ... but the tramTMB object is set up with using them
chkeq(ip$beta, mod$tmb_obj$env$par_checked, check.attributes = FALSE)

## -- Data can come from the gobal environment
data("sleepstudy", package = "lme4")
fit_lm1 <- LmME(Reaction ~ Days + (Days || Subject), data = sleepstudy)
attach(sleepstudy)
fit_lm2 <- LmME(Reaction ~ Days + (Days || Subject))
chkeq(logLik(fit_lm1), logLik(fit_lm2))

options(oldopt)
