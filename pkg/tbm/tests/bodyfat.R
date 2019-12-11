
library("tbm")
library("tram")
library("partykit")

set.seed(29)

data("bodyfat", package = "TH.data")

mf <- as.mlt(Colr(DEXfat ~ 1, data = bodyfat, order = 5))
logLik(mf)

Mstop <- 50

fd <- cv(rep(1, NROW(bodyfat)), type = "kfold", B = 2)

bctrl <- boost_control(nu = .1, trace = FALSE, mstop = Mstop)

tctrl <- ctree_control(minsplit = 8, minbucket = 5, mincriterion = 0,
                       maxdepth = 4, splittest = TRUE, 
                       testtype = "Teststatistic")

bf_t <- ctmboost(model = mf, formula = DEXfat ~ ., data = bodyfat, 
                 method = quote(mboost::blackboost), control = bctrl, 
                 tree_control = tctrl)
logLik(bf_t)

bf_ctm <- ctmboost(model = mf, formula = DEXfat ~ ., data = bodyfat, 
                   control = bctrl)
logLik(bf_ctm)
table(selected(bf_ctm))

bf_dr <- ctmboost(model = mf, formula = DEXfat ~ ., data = bodyfat,
                  baselearner = "bols", control = bctrl)
logLik(bf_dr)
table(selected(bf_dr))

bf_st <- stmboost(model = mf, formula = DEXfat ~ ., data = bodyfat, 
                  method = quote(mboost::blackboost), tree_control = tctrl)
round(logLik(bf_st), 1)

bf_shift <- stmboost(model = mf, formula = DEXfat ~ ., data = bodyfat, 
                     method = quote(mboost::gamboost))
logLik(bf_shift)
table(selected(bf_shift))

bf_lin <- stmboost(model = mf, formula = DEXfat ~ . - 1, data = bodyfat, 
                   method = quote(mboost:::glmboost.formula))
logLik(bf_lin)
table(selected(bf_lin))

mf2 <- Lm(DEXfat ~ 1, data = bodyfat)

bf_lin2 <- ctmboost(model = mf2, formula = DEXfat ~ ., data = bodyfat)
logLik(bf_lin2$model, parm = coef(bf_lin2))
table(selected(bf_lin2))

### test against L_2 glmboost
m <- Lm(DEXfat ~ 1, data = bodyfat, fixed = c("DEXfat" = 1))
bf_1 <- stmboost(model = m, formula = DEXfat ~ 0 + ., data = bodyfat, 
                  control = bctrl,
                  method = quote(mboost:::glmboost.formula), 
                  mltargs = list(fixed = c("DEXfat" = 1)))
bf_2 <- glmboost(DEXfat ~ ., data = bodyfat, offset = mean(bodyfat$DEXfat),
                  control = bctrl)
stopifnot(max(abs(mboost:::coef.glmboost(bf_1) - coef(bf_2)[-1])) < 
          sqrt(.Machine$double.eps))
r <- risk(bf_1)
stopifnot(r[length(r)] + logLik(bf_1) < sqrt(.Machine$double.eps))

stopifnot(max(abs(-nuisance(bf_1) + mboost:::predict.glmboost(bf_1) - predict(bf_2))) < 
          sqrt(.Machine$double.eps))

