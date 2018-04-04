
library("tram")
library("survival")
library("tbm")
library("lattice")

data("GBSG2", package = "TH.data")
GBSG2y <- numeric_var("y", support = c(100.0, max(GBSG2$time)), 
                      bounds = c(0, Inf))
GBSG2$y <- with(GBSG2, Surv(time, cens))
B_GBSG2y <- Bernstein_basis(var = GBSG2y, order = 6, ui = "increasing")
var_a <- numeric_var("age", support = range(GBSG2$age))
B_age <- Bernstein_basis(var_a, order = 3)
b_horTh <- as.basis(GBSG2$horTh)
ctm_GBSG2_8 <- ctm(B_GBSG2y, 
                   interacting = b(horTh = b_horTh, age = B_age), 
                   todistr = "MinExtrVal")
mlt_GBSG2_8  <- mlt(ctm_GBSG2_8, data = GBSG2)
logLik(mlt_GBSG2_8)

ctm_GBSG2_8_2 <- ctm(B_GBSG2y, 
                   interacting = B_age, 
                   todistr = "MinExtrVal")
mlt_GBSG2_8_2  <- mlt(ctm_GBSG2_8_2, data = GBSG2)
logLik(mlt_GBSG2_8_2)


s <- mkgrid(mlt_GBSG2_8, 100)
nd <- expand.grid(s)
nd$s <- c(predict(mlt_GBSG2_8_2, newdata = s, type = "survivor"))
contourplot(s ~ age + y | horTh, data = nd, at = 1:9 / 10,
            ylab = "Survival time (days)", xlab = "Age (years)",
            scales = list(x = list(alternating = c(1, 1))))

nd$s <- c(predict(mlt_GBSG2_8, newdata = s, type = "survivor"))
contourplot(s ~ age + y | horTh, data = nd, at = 1:9 / 10,
            ylab = "Survival time (days)", xlab = "Age (years)",
            scales = list(x = list(alternating = c(1, 1))))



m <- Coxph(y ~ 1, data = GBSG2)

tb <- ctmboost(m, y ~ bbs(age),# + bbs(age, by = horTh), 
          control = boost_control(nu = .05, trace = TRUE), 
          data = GBSG2)[500]

fd <- cv(model.weights(tb), type = "sub", B = 5, strata = GBSG2$horTh)
#cv <- cvrisk(tb, folds = fd)
#plot(cv)
#tb <- tb[mstop(cv)]

logLik(m, parm = coef(tb))

s$y <- NULL
nd2 <- expand.grid(s)

nd$s2 <- c(predict(tb, newdata = nd2, type = "survivor", K = 100))
contourplot(s2 ~ age + y | horTh, data = nd, at = 1:9 / 10,
            ylab = "Survival time (days)", xlab = "Age (years)",
            scales = list(x = list(alternating = c(1, 1))))


tb <- ctmboost(m, y ~ bbs(age) + bbs(age, by = horTh), 
          control = boost_control(nu = .01, trace = TRUE), 
          data = GBSG2)[1000]

cv <- cvrisk(tb, folds = fd)
plot(cv)
tb <- tb[mstop(cv)]

logLik(m, parm = coef(tb))

s$y <- NULL
nd2 <- expand.grid(s)

nd$s2 <- c(predict(tb, newdata = nd2, type = "survivor", K = 100))
contourplot(s2 ~ age + y | horTh, data = nd, at = 1:9 / 10,
            ylab = "Survival time (days)", xlab = "Age (years)",
            scales = list(x = list(alternating = c(1, 1))))


tc <- partykit:::ctree_control(maxdepth = 3, mincriterion = 0)

tb2 <- ctmboost(m, y ~ horTh + age + menostat + tsize + tgrade + pnodes + 
              progrec + estrec, data = GBSG2,
              control = boost_control(nu = .01), method = quote(mboost::blackboost),
              tree_control = tc)[1000]

logLik(m, parm = coef(tb2))

nd <- GBSG2
nd$horTh[] <- "yes"
tyes <- predict(tb2, newdata = nd, type = "trafo")
nd$horTh[] <- "no"
tno <- predict(tb2, newdata = nd, type = "trafo")

matplot(tyes - tno, type = "l")

summary(Coxph(y ~ horTh + age + menostat + tsize + tgrade + pnodes + 
              progrec + estrec, data = GBSG2))

