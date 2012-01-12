
library("ctm")
library("quantreg")

if (FALSE) {

### sorry, these data are not publically available
load("fetus.Rda")

### linear model
lmod <- lm(birthweight ~ volabdo + hccalc +
	   volos +fe + bip + I(bip^2), data = fetus)

### predictions
pred1 <- as.data.frame(predict(lmod, interval = "prediction"))
pred1$birthweight <- fetus$birthweight
pred1 <- pred1[order(pred1$fit),]
names(pred1) <- c("q50", "q10", "q90", "birthweight")
pred1 <- pred1[, c("q10", "q50", "q90", "birthweight")]

### linear quantile regression
xvar <- colnames(fetus)[-1]
fm <- paste("birthweight ~ ", paste(xvar, collapse = "+"))
rq10 <- rq(as.formula(fm), data = fetus, tau = .1)
rq50 <- rq(as.formula(fm), data = fetus, tau = .5)
rq90 <- rq(as.formula(fm), data = fetus, tau = .9)

pred3 <- data.frame(q10 = fitted(rq10),
                 q50 = fitted(rq50),
                 q90 = fitted(rq90))
pred3$birthweight <- fetus$birthweight  
pred3 <- pred3[order(pred3$q50),]

### additive quantile regression
fm <- paste("birthweight ~ ", paste("bbs(", xvar, ", df = 5)", collapse = "+"))
ctrl <- boost_control(nu = .9)
rqss10 <- mboost(as.formula(fm), data = fetus, family = QuantReg(tau = .1), control = ctrl)[2000]
rqss50 <- mboost(as.formula(fm), data = fetus, family = QuantReg(tau = .5), control = ctrl)[2000]
rqss90 <- mboost(as.formula(fm), data = fetus, family = QuantReg(tau = .9), control = ctrl)[2000]
cv10 <- cvrisk(rqss10)
cv50 <- cvrisk(rqss50)
cv90 <- cvrisk(rqss90)

pred4 <- data.frame(q10 = fitted(rqss10),
                    q50 = fitted(rqss50),
                    q90 = fitted(rqss90))
pred4$birthweight <- fetus$birthweight  
pred4 <- pred4[order(pred4$q50),]

### ctm
fm <- paste("bbs(birthweight, df = 2.05) ~ ", paste("bbs(", xvar, ", df = 2.05)", 
            collapse = "+"))
dmod <- ctm(as.formula(fm), data = fetus, family = Binomial(link = "probit"),
               control = boost_control(nu = 0.2, mstop = 500))
# tune(dmod, alpha = .0001)
(cv <- cvrisk(dmod))
dmod[mstop(cv)]

p <- predict(dmod, newdata = fetus, anno = TRUE, type = "response")
pred2 <- tapply(1:NROW(p), p$ID, function(i) qest(p[i,]))
pred2 <- matrix(unlist(pred2), ncol = 3, byrow = TRUE)
pred2 <- data.frame(pred2)
colnames(pred2) <- c("q10", "q50", "q90")
pred2$birthweight <- fetus$birthweight
pred2 <- pred2[order(pred2$q50),]

### tune and predict won't work because variable ONE is missing
if (FALSE) {
fetus$EINS <- 1
fm <- paste("bols(birthweight, intercept = FALSE, df = 1) + 
             bols(ONE, intercept = FALSE, df = 1) + 
             bbs(birthweight, center = TRUE, df = 1) ~ ", 
paste(c("bols(EINS, intercept = FALSE, df = 1)", 
      paste("bols(", xvar, ", intercept = FALSE, df = 1)"),
      paste("bbs(", xvar, ", center = TRUE, df = 1)")), collapse = "+"))
fm <- as.formula(fm)
dmod2 <- ctm(fm, data = fetus, family = Binomial(link = "probit"))
tune(dmod2, alpha = .0001)

p <- predict(dmod2, newdata = fetus, anno = TRUE, type = "response")
pred5 <- tapply(1:NROW(p), p$ID, function(i) qest(p[i,]))
pred5 <- matrix(unlist(pred5), ncol = 3, byrow = TRUE)
pred5 <- data.frame(pred5)
colnames(pred5) <- c("q10", "q50", "q90")
pred5$birthweight <- fetus$birthweight
pred5 <- pred5[order(pred5$q50),]

}


pred1$model <- "lm"
pred1$indx <- 1:nrow(pred1)
pred2$model <- "ctm"
pred2$indx <- 1:nrow(pred2)
pred3$model <- "rq"
pred3$indx <- 1:nrow(pred3)
pred4$model <- "rqss"
pred4$indx <- 1:nrow(pred4)

pred <- rbind(pred1, pred2, pred3, pred4)
pred$model <- factor(pred$model)

}

save(pred, file = "ex_fetus.Rda")

