
library("ctm")

load("DVC.Rda")

DVC <- data.frame(number = c(rowSums(x06), rowSums(x09)),
                  days = c(1:365, 1:365),
                  year = factor(c(rep("2006", 365), rep("2009", 365))))

mod <- ctm(bbs(number, df = 3, differences = 1) ~ 
                bbs(days, df = 3, differences = 1, cyclic = TRUE)
              + bbs(days, by = year, differences = 1, df = 3, cyclic = TRUE),
              data = DVC, family = Binomial(link = "probit"),
              control = boost_control(nu = 0.5, mstop = 250))

folds <- cv(model.weights(mod), strata = DVC$year)
(cv <- cvrisk(mod, folds = folds))
mod[mstop(cv)]

p <- predict(mod, newdata = DVC, annotate = TRUE, type = "response")
DVC$ID <- 1:nrow(DVC)
DVC$days <- c(seq(as.Date("2006-01-01"), length.out=365, by="1 day"),
              seq(as.Date("2009-01-01"), length.out=365, by="1 day"))
p <- merge(p, DVC, by = "ID")

save(p, file = "ex_DVC.Rda")
