
library("ctm")   
library("np")

set.seed(290875)

data("fattyacid", package = "multcomp")
xdf <- subset(fattyacid, complete.cases(fattyacid))

fit <- c()
grid <- seq(from = min(xdf$FA), to = max(xdf$FA), length = 100)
y <- rep(grid, rep(nlevels(xdf$PE), length(grid)))
x <- rep(unique(xdf$PE), length(grid))
ctrl <- boost_control(mstop = 250)

fa <- lapply(xdf$FA, function(y)
    with(xdf, data.frame(FAtmp = FA > y, PE = PE, FA = y)))
fa <- do.call("rbind", fa)
fa$FAtmp <- factor(fa$FAtmp)
m <- glm(FAtmp ~ PE + PE:FA - 1, data = fa, family = binomial(link = "probit"))
nfa <- expand.grid(PE = unique(xdf$PE), FA = grid)
p <- predict(m, newdata = nfa, type = "response")

fit <- rbind(fit, data.frame(p = p, y = y, x = x, what = "glm",
                             stringsAsFactors = FALSE))

### df = 6 groups x 3 df each
mod <- ctm(bbs(FA, df = 18) ~ bols(PE, intercept = FALSE, 
	                                lambda = 0), data = xdf,
              family = Binomial(link = "probit"), monotone = FALSE,
              control = ctrl)
# tune(mod)
folds <- cv(model.weights(mod), strata = xdf$PE)
(cv <- cvrisk(mod, folds = folds))
mod[mstop(cv)]
p <- predict(mod, newdata = data.frame(PE = unique(xdf$PE)), 
             y = grid, type = "response")
fit <- rbind(fit, data.frame(p = p, y = y, x = x, what = "bbs probit", 
                             stringsAsFactors = FALSE))


mod <- ctm(bols(FA, df = 8) ~ bols(PE, intercept = FALSE, 
                                      lambda = 0), data = xdf,
              family = Binomial(link = "probit"), monotone = FALSE,
              control = ctrl)
folds <- cv(model.weights(mod), strata = xdf$PE)
(cv <- cvrisk(mod, folds = folds))
mod[mstop(cv)]
# tune(mod)
p <- predict(mod, newdata = data.frame(PE = unique(xdf$PE)), 
             y = grid, type = "response")
fit <- rbind(fit, data.frame(p = p, y = y, x = x, what = "bols probit", 
                             stringsAsFactors = FALSE))


mod <- npcdist(FA ~ PE, data = xdf)
p <- predict(mod, newdata = data.frame(FA = y, PE = x))
fit <- rbind(fit, data.frame(p = p, y = y, x = x, what = "npcdist", 
                             stringsAsFactors = FALSE))

save.image(file = "ex_fattyacid.Rda")
