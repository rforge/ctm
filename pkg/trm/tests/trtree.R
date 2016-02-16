
library("trm")
library("partykit")
library("survival")
data("GBSG2", package = "TH.data")

yvar <- numeric_var("y", support = c(100, 2000), bounds = c(0, Inf))
By <- Bernstein_basis(yvar, order = 5, ui = "incre")
m <- ctm(response = By, todistr = "MinExt")
GBSG2$y <- with(GBSG2, Surv(time, cens))
mlmod <- mlt(m, data = GBSG2, scale = TRUE, check = FALSE, 
             checkGrad = FALSE, gtol = 1e-3)

ctrl <- ctree_control(minsplit = 40)
tt <- trtree(mlmod, ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2, control = ctrl)

tt2 <- trtree(mlmod, ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2, control = ctrl, 
    modelsplit = TRUE)

plot(tt, tp_args = list(type = "surv"))
plot(tt2, tp_args = list(type = "surv"))

ct <- ctree(Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2, control = ctrl) 

plot(ct)

