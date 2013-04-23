
source("simfuns.R")

set.seed(290875)

pnon <- 0
xdf <- dgp(c(100, 100), pnon = pnon)
test <- dgp(pnon = pnon)

boost <- Fboost(xdf)
pt <- pdf(truth, test)
pb <- pdf(boost, test)

q <- function(tau = 0.5) {

    mt <- apply(pt, 2, function(x) {
        af <- approxfun(ys, x)
        f <- function(y) (af(y) - tau)^2
        optimize(f, range(ys))$minimum
    })

    mb <- apply(pb, 2, function(x) {
        af <- approxfun(ys, x)
        f <- function(y) (af(y) - tau)^2
        optimize(f, range(ys))$minimum
    })

    qboost <- mboost(y ~ bbs(x1, by = x2, df = 4) + bbs(x1, df = 4), 
                     data = xdf, family = QuantReg(tau = tau))
    qboost[500]
    cv <- cvrisk(qboost, folds = matrix(as.integer(xdf$ng == "learn")))
    qboost[mstop(cv)]
    mq <- predict(qboost, newdata = test)

    ret <- data.frame(mt = mt, mb = mb, mq = mq, tau = tau)
}

q5 <- q(0.5)
q75 <- q(0.75)
q9 <- q(0.9)

q <- rbind(q5, q75, q9)
q$tau <- as.factor(paste(q$tau, "quantile"))
save(q, file = "qret.Rda")
