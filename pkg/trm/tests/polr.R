
library("trm")
library("multcomp")
library("MASS")

options(contrasts = c("contr.treatment", "contr.poly"))
fm <- Sat ~ Infl + Type
house.plr <- polr(fm, weights = Freq, data = housing)

house.sltm <- sltm(fm, weights = Freq, data = housing)

summary(house.plr, digits = 3)
summary(house.sltm, digits = 3)  

summary(update(house.plr, method = "probit", Hess = TRUE), digits = 3)
summary(sltm(fm, weights = Freq, data = housing, method = "probit"))

summary(update(house.plr, method = "cloglog", Hess = TRUE), digits = 3)
summary(sltm(fm, weights = Freq, data = housing, method = "cloglog"))

max(abs(predict(house.plr, housing, type = "p") - 
        t(predict(house.sltm, type = "dens"))))

pr <- profile(house.plr)
confint(pr)

plot(house.sltm)
