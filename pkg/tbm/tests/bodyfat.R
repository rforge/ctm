
library("tbm")
library("tram")

set.seed(29)

layout(matrix(1:4, ncol = 4))

data("bodyfat", package = "TH.data")

mf <- as.mlt(Colr(DEXfat ~ 1, data = bodyfat, order = 5))
logLik(mf)

fd <- cv(rep(1, NROW(bodyfat)), type = "kfold", B = 5)

bf_ctm <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bbs")
ms <- cvrisk(bf_ctm, folds = fd)
plot(ms)
bf_ctm <- bf_ctm[mstop(ms)]
logLik(mf, parm = coef(bf_ctm))
table(selected(bf_ctm))

bf_dr <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bols")
ms <- cvrisk(bf_dr, folds = fd)
plot(ms)
bf_dr <- bf_dr[mstop(ms)]
logLik(mf, parm = coef(bf_dr))
table(selected(bf_dr))

bf_shift <- tbm(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bbs", 
                gradient = "shift")
ms <- cvrisk(bf_shift, folds = fd)
plot(ms)
bf_shift <- bf_shift[mstop(ms)]
logLik(bf_shift$model, parm = coef(bf_shift))
table(selected(bf_shift))

bf_lin <- tbm(model = mf, formula = DEXfat ~ . - 1, data = bodyfat, baselearner = "bols", 
                gradient = "shift")
ms <- cvrisk(bf_lin, folds = fd)
plot(ms)
bf_lin <- bf_lin[mstop(ms)]
logLik(bf_lin$model, parm = coef(bf_lin))
table(selected(bf_lin))


col <- rgb(.1, .1, .1, .1)
matplot(predict(bf_ctm, type = "distribution"), type = "l", main = "CTM", col =
col, lty = 1)
matplot(predict(bf_dr, type = "distribution"), type = "l", main = "Distribution Regression", col = col, lty = 1)
matplot(predict(bf_shift, type = "distribution"), type = "l", main = "GAM", col = col, lty = 1)
matplot(predict(bf_lin, type = "distribution"), type = "l", main = "TRAM", col = col, lty = 1)


