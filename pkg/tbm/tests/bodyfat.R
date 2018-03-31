
library("tbm")
library("tram")

set.seed(29)

layout(matrix(1:6, ncol = 3, byrow = TRUE))

data("bodyfat", package = "TH.data")

mf <- as.mlt(Colr(DEXfat ~ 1, data = bodyfat, order = 5))
logLik(mf)

Mstop <- 110

fd <- cv(rep(1, NROW(bodyfat)), type = "kfold", B = 5)

schwarzboost <- tbm:::schwarzboost

bf_t <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bbaum")[Mstop]
ms <- cvrisk(bf_t, folds = fd)
plot(ms, main = "CTM-Baum")
bf_t <- bf_t[mstop(ms)]
logLik(mf, parm = coef(bf_t, which = 1))


bf_ctm <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bbs")[Mstop]
ms <- cvrisk(bf_ctm, folds = fd)
plot(ms, main = "CTM-Add")
bf_ctm <- bf_ctm[mstop(ms)]
logLik(mf, parm = coef(bf_ctm))
table(selected(bf_ctm))

bf_dr <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bols")[Mstop]
ms <- cvrisk(bf_dr, folds = fd)
plot(ms, main = "Distr Reg")
bf_dr <- bf_dr[mstop(ms)]
logLik(mf, parm = coef(bf_dr))
table(selected(bf_dr))

bf_st <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bbaum", 
                gradient = "shift")[Mstop]
ms <- cvrisk(bf_st, folds = fd)
plot(ms, main = "Baum")
bf_st <- bf_st[mstop(ms)]
logLik(bf_st$model, parm = coef(bf_st))

bf_shift <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bbs", 
                gradient = "shift")[Mstop]
ms <- cvrisk(bf_shift, folds = fd)
plot(ms, main = "GAM")
bf_shift <- bf_shift[mstop(ms)]
logLik(bf_shift$model, parm = coef(bf_shift))
table(selected(bf_shift))

bf_lin <- tbm(model = mf, formula = DEXfat ~ . - 1, data = bodyfat, baselearner = "bols", 
                gradient = "shift")[Mstop]
ms <- cvrisk(bf_lin, folds = fd)
plot(ms, main = "TRAM")
bf_lin <- bf_lin[mstop(ms)]
logLik(bf_lin$model, parm = coef(bf_lin))
table(selected(bf_lin))


col <- rgb(.1, .1, .1, .1)
matplot(predict(bf_t, type = "distribution"), type = "l", main = "CTM-Baum", col = col, lty = 1)
matplot(predict(bf_ctm, type = "distribution"), type = "l", main = "CTM-Add", col = col, lty = 1)
matplot(predict(bf_dr, type = "distribution"), type = "l", main = "Distribution Regression", col = col, lty = 1)
matplot(predict(bf_st, newdata = bodyfat, type = "distribution"), type = "l", main = "Baum", col = col, lty = 1)
matplot(predict(bf_shift, newdata = bodyfat, type = "distribution"), type = "l", main = "GAM", col = col, lty = 1)
matplot(predict(bf_lin, newdata = bodyfat, type = "distribution"), type = "l", main = "TRAM", col = col, lty = 1)
