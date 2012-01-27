
library("ctm")   
library("np")
set.seed(290875)

data("Italy", package = "np")
Italy$year <- with(Italy, as.numeric(as.character(year)))
mod <- ctm(bbs(gdp, df = 3) ~ 
              bbs(year, df = 3), data = Italy, 
              family = Binomial(link = "probit"),
              control = boost_control(trace = TRUE, nu = .2, mstop = 500),
              ngrid = 50)

(cv <- cvrisk(mod))
mod[mstop(cv)]

pr <- expand.grid(year = min(Italy$year):max(Italy$year),
                  gdp = with(Italy, seq(from = min(gdp), to = max(gdp),
                                        length = 25)))
pr$prob <- predict(mod, newdata = pr, y = NULL, type = "response")
pr$model <- "ctm"

npp <- pr
nd <- npcdist(gdp ~ year, data = Italy)
npp$model <- "np"
npp$prob <- predict(nd, newdata = npp)
pr <- rbind(pr, npp)

save(pr, file = "ex_italy.Rda")
