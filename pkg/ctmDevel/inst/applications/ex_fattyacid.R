
source("setup.R")
library("np")

set.seed(290875)

data("fattyacid", package = "multcomp")
xdf <- subset(fattyacid, complete.cases(fattyacid))

fit <- c()
grid <- seq(from = min(xdf$FA), to = max(xdf$FA), length = 100)
y <- rep(grid, rep(nlevels(xdf$PE), length(grid)))
x <- rep(unique(xdf$PE), length(grid))
ctrl <- boost_control(mstop = 250, trace = TRUE)

fa <- lapply(xdf$FA, function(y)
    with(xdf, data.frame(FAtmp = FA > y, PE = PE, FA = y)))
fa <- do.call("rbind", fa)
fa$FAtmp <- factor(fa$FAtmp)
m <- glm(FAtmp ~ PE + PE:FA - 1, data = fa, family = binomial(link = "probit"))
nfa <- expand.grid(PE = unique(xdf$PE), FA = grid)
p <- predict(m, newdata = nfa, type = "response")

fit <- rbind(fit, data.frame(p = p, y = y, x = x, what = "glm",
                             stringsAsFactors = FALSE))

xdf <- cbind(xdf, model.matrix(~ PE - 1, data = xdf))

### df = 6 groups x 3 df each
### use separate base-learners for each level -> separate smoothness
### we see slight violations of monotonicity outside the range of
### FA for some groups
mod <- mctm(bbs(FA, df = 3) ~ 
    bols(PEPE3, intercept = FALSE, df = 1) + 
    bols(PEPE4, intercept = FALSE, df = 1) + 
    bols(PEPE5, intercept = FALSE, df = 1) + 
    bols(PEPE6, intercept = FALSE, df = 1) + 
    bols(PEPE7, intercept = FALSE, df = 1) + 
    bols(PEPE9, intercept = FALSE, df = 1), 
                                    data = xdf, ngrid = sort(unique(xdf$FA)),
              family = Binomial(link = "probit"), 
              control = boost_control(nu = .4, mstop = 1000, trace = TRUE))
# mod[1000]
# tune(mod)
# folds <- cv(model.weights(mod), strata = xdf$PE)
#(cv <- cvrisk(mod, folds = folds))
#mod[mstop(cv)]
I <- diag(nlevels(xdf$PE))
colnames(I) <- paste("PE", as.character(unique(xdf$PE)), sep = "")
I <- as.data.frame(I)
p <- predict(mod, newdata = I,
             y = grid, type = "response")
fit <- rbind(fit, data.frame(p = p, y = y, x = x, what = "bbs probit",
                             stringsAsFactors = FALSE))


mod <- ctm(bols(FA, df = 2) ~ bols(PE, df = 6, contrasts.arg = "contr.dummy"),
              data = xdf,
              family = Binomial(link = "probit"), monotone = FALSE,
              control = ctrl)
folds <- cv(model.weights(mod), strata = xdf$PE)
#(cv <- cvrisk(mod, folds = folds))
#mod[mstop(cv)]
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
