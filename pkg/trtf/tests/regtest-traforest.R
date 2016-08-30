
library("trtf")
library("survival")
data("GBSG2", package = "TH.data")

yvar <- numeric_var("y", support = c(100, 2000), bounds = c(0, Inf))
By <- Bernstein_basis(yvar, order = 5, ui = "incre")
m <- ctm(response = By, todistr = "MinExt")
GBSG2$y <- with(GBSG2, Surv(time, cens))

tf <- traforest(m, formula = y ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2, 
    control = ctree_control(splitstat = "quad", teststat = "quad",
                    testtype = "Teststatistic", mincriterion = 1, minbucket = 50), 
    ntree = 50, trace = TRUE, cores = 4)

w <- predict(tf, newdata = GBSG2[1:3,], type = "weights")

cf <- coef(mlt(m, data = GBSG2))
coef(m1 <- mlt(m, data = GBSG2, weights = w[,1], theta = cf))
coef(m2 <- mlt(m, data = GBSG2, weights = w[,2], theta = cf))
coef(m3 <- mlt(m, data = GBSG2, weights = w[,3], theta = cf))

layout(matrix(1:3, nr = 1))
plot(m1, newdata = data.frame(1), type = "survivor")
plot(m2, newdata = data.frame(1), type = "survivor", add = TRUE)
plot(m3, newdata = data.frame(1), type = "survivor", add = TRUE)

cmod <- coxph(Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2)

plot(survfit(cmod, newdata = GBSG2[1:3,]))

sf <-  cforest(formula = y ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2, 
    control = ctree_control(splitstat = "quad", teststat = "quad",
                    testtype = "Teststatistic", mincriterion = 1, minbucket = 50),
    ntree = 50, trace = TRUE)

w <- predict(sf, newdata = GBSG2[1:3,], type = "weights")

cf <- coef(mlt(m, data = GBSG2))
coef(m1 <- mlt(m, data = GBSG2, weights = w[,1], theta = cf))
coef(m2 <- mlt(m, data = GBSG2, weights = w[,2], theta = cf))
coef(m3 <- mlt(m, data = GBSG2, weights = w[,3], theta = cf))

plot(m1, newdata = data.frame(1), type = "survivor")
plot(m2, newdata = data.frame(1), type = "survivor", add = TRUE)
plot(m3, newdata = data.frame(1), type = "survivor", add = TRUE)
