# Test for cvl_tramnet

library("survival")
library("tramnet")
library("TH.data")
options(digits = 3)

set.seed(241068)
data("GBSG2", package = "TH.data")
X <- 1 * matrix(GBSG2$horTh == "yes", ncol = 1)
colnames(X) <- "horThyes"
GBSG2$surv <- with(GBSG2, Surv(time, cens))
m <- Coxph(surv ~ 1, data = GBSG2)
mt <- tramnet(model = m, x = X, lambda = 0, alpha = 0)
mc <- Coxph(surv ~ horTh, data = GBSG2)
cvl_tramnet(mt, fold = 2, lambda = c(0, 1), alpha = c(0, 1))
