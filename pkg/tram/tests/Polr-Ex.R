
library("MASS")
library("tram")

(house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data =
housing))

summary(house.plr)

(h <- Polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing))

summary(h)




