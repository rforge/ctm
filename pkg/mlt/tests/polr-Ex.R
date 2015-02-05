
library("MASS")

mp <- polr(Sat ~ Infl, weights = Freq, data = housing)

library("mlt")

s <- as.basis(~ Infl, data = housing, remove_intercept = TRUE)
r <- as.basis(~ Sat, data = housing, remove_intercept = TRUE,
              contrasts.arg = list(Sat = function(n) 
                  contr.treatment(n, base = 3)),
              ui = diff(Diagonal(2)), ci = 0)

m <- model(r, shift = s)

mod <- mlt(m, housing, weights = housing$Freq, todist = "Logi")

logLik(mp)
logLik(mod)

coef(mp)
mp$zeta
coef(mod)
