
library("tbm")
library("tram")

set.seed(29)

layout(matrix(1:6, ncol = 3, byrow = TRUE))

data("bodyfat", package = "TH.data")

mf <- as.mlt(Colr(DEXfat ~ 1, data = bodyfat, order = 5))
logLik(mf)

Mstop <- 500

fd <- cv(rep(1, NROW(bodyfat)), type = "kfold", B = 2)

bf_t <- ctmboost(model = mf, formula = DEXfat ~ ., data = bodyfat, method = quote(mboost::blackboost))[Mstop]
ms <- cvrisk(bf_t, folds = fd)
plot(ms, main = "CTM-Baum")
bf_t <- bf_t[mstop(ms)]
logLik(mf, parm = coef(bf_t))

bf_ctm <- ctmboost(model = mf, formula = DEXfat ~ ., data = bodyfat)[Mstop]
ms <- cvrisk(bf_ctm, folds = fd)
plot(ms, main = "CTM-Add")
bf_ctm <- bf_ctm[mstop(ms)]
logLik(mf, parm = coef(bf_ctm))
table(selected(bf_ctm))

bf_dr <- ctmboost(model = mf, formula = DEXfat ~ ., data = bodyfat, baselearner = "bols")[Mstop]
ms <- cvrisk(bf_dr, folds = fd)
plot(ms, main = "Distr Reg")
bf_dr <- bf_dr[mstop(ms)]
logLik(mf, parm = coef(bf_dr))
table(selected(bf_dr))

bf_st <- tramboost(model = mf, formula = DEXfat ~ ., data = bodyfat, method =
quote(mboost::blackboost))[Mstop]
ms <- cvrisk(bf_st, folds = fd)
plot(ms, main = "Baum")
bf_st <- bf_st[mstop(ms)]
logLik(bf_st$model, parm = coef(bf_st))

bf_shift <- tramboost(model = mf, formula = DEXfat ~ ., data = bodyfat, method = quote(mboost::gamboost))[Mstop]
ms <- cvrisk(bf_shift, folds = fd)
plot(ms, main = "GAM")
bf_shift <- bf_shift[mstop(ms)]
logLik(bf_shift$model, parm = coef(bf_shift))
table(selected(bf_shift))

bf_lin <- tramboost(model = mf, formula = DEXfat ~ . - 1, data = bodyfat, method = quote(mboost:::glmboost.formula))[Mstop]
ms <- cvrisk(bf_lin, folds = fd)
plot(ms, main = "TRAM")
bf_lin <- bf_lin[mstop(ms)]
logLik(bf_lin$model, parm = coef(bf_lin))
table(selected(bf_lin))


mf2 <- Lm(DEXfat ~ 1, data = bodyfat)

bf_lin2 <- ctmboost(model = mf2, formula = DEXfat ~ ., data = bodyfat)[Mstop]
ms <- cvrisk(bf_lin2, folds = fd)
plot(ms, main = "TRAM-exp")
bf_lin2 <- bf_lin2[mstop(ms)]
logLik(bf_lin2$model, parm = coef(bf_lin2))
table(selected(bf_lin2))

X11()
layout(matrix(1:6, ncol = 3, byrow = TRUE))

col <- rgb(.1, .1, .1, .1)
matplot(predict(bf_t, newdata = bodyfat, type = "distribution"), type = "l", main = "CTM-Baum", col = col, lty = 1)
matplot(predict(bf_ctm, newdata = bodyfat, type = "distribution"), type = "l", main = "CTM-Add", col = col, lty = 1)
matplot(predict(bf_dr, newdata = bodyfat, type = "distribution"), type = "l", main = "Distribution Regression", col = col, lty = 1)
matplot(predict(bf_st, newdata = bodyfat, type = "distribution"), type = "l", main = "Baum", col = col, lty = 1)
matplot(predict(bf_shift, newdata = bodyfat, type = "distribution"), type = "l", main = "GAM", col = col, lty = 1)
matplot(predict(bf_lin, newdata = bodyfat, type = "distribution"), type = "l", main = "TRAM", col = col, lty = 1)



