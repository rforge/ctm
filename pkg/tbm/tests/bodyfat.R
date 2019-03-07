
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

tctrl <- ctree_control(minsplit = 2, minbucket = 1, mincriterion = 0,
                       maxdepth = 5, splittest = TRUE, 
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
logLik(bf_st)

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
