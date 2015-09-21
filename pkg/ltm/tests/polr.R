
library("ltm")
library("multcomp")
library("MASS")

options(contrasts = c("contr.treatment", "contr.poly"))
fm <- Sat ~ Infl + Type
house.plr <- polr(fm, weights = Freq, data = housing)

house.ltm <- ltm(fm, weights = Freq, data = housing)

summary(house.plr, digits = 3)
summary(house.ltm, digits = 3)  

summary(update(house.plr, method = "probit", Hess = TRUE), digits = 3)
summary(ltm(fm, weights = Freq, data = housing, method = "probit"))

summary(update(house.plr, method = "cloglog", Hess = TRUE), digits = 3)
summary(ltm(fm, weights = Freq, data = housing, method = "cloglog"))

max(abs(predict(house.plr, housing, type = "p") - 
        t(predict(house.ltm, type = "dens"))))

pr <- profile(house.plr)
confint(pr)

plot(house.ltm)
