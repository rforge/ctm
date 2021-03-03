## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4, warn = -1)
set.seed(100)

library("tramME")
library("survival")

if (FALSE) { ## NOTE: do not run it for now, because it requires a non-CRAN version of tram

## -- compare to tram
m1 <- Aareg(Surv(tstart, tstop, status) | age + treat ~ 1, data = cgd, order = 6)
stopifnot(m1$convergence == 0)
m2 <- AaregME(Surv(tstart, tstop, status) | age + treat ~ 1, data = cgd, order = 6)
stopifnot(m2$opt$convergence == 0)
chkeq(logLik(m1), logLik(m2), check.attributes = FALSE)
chkeq(coef(m1, with_baseline = TRUE), coef(m2, with_baseline = TRUE), tol = 1e-4)

## -- mixed-effects version
if (.run_test) {
  m3 <- AaregME(Surv(tstart, tstop, status) | age + treat + sex ~ 1 + (1 | id),
                data = cgd, order = 8)
  stopifnot(m3$opt$convergence == 0)
}

}
