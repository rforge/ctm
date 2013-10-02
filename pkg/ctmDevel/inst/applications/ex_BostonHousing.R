
source("setup.R")
set.seed(290875)

library("spdep")
library("Matrix")

data("boston")

nb <- nb2mat(boston.soi, zero.policy=TRUE, style="B")
nb <- -nb
diag(nb) <- -rowSums(nb)
colnames(nb) <- rownames(nb) <- 1:nrow(boston.c)
# nb <- Matrix(nb)

weights <- as.double(with(boston.c, CMEDV < max(CMEDV)))

nm <- names(boston.c)[c(8, 9, 10, 12:20)]
for (n in nm) boston.c[[n]] <- as.vector(scale(boston.c[[n]]))

bk <- range(boston.c$CMEDV) - sd(boston.c$CMEDV) * c(1, -1)

boston.c$EINS <- boston.c$ONE <- 1.0

BH <- boston.c[c(nm, "ONE", "CMEDV", "EINS")]
BH$ID <- factor(1:nrow(BH))

if (FALSE) {
fm <- paste("bols(CMEDV, intercept = FALSE, df = 1) + 
###             bols(ONE, intercept = FALSE, df = 1) + 
             bbs(CMEDV, df = 1, boundary.knots = bk, center = TRUE) ~ ", 
             paste(c("bols(EINS, intercept = FALSE, df = 1)", 
                    paste("bols(", nm, ", df = 1, intercept = FALSE)"), 
                    paste("bbs(", nm, ", df = 1, center = TRUE)")), 
            collapse = "+"))
fm <- as.formula(fm)

options(mboost_eps = 10e-6)

mod <- ctm(fm, data = BH, 
              family = Binomial(link = "probit"),  monotone = FALSE,
              control = boost_control(mstop = 50, trace = TRUE), 
              ngrid = 25, weights = weights, 
              constant = "bmrf(ID, bnd = nb, df = 1)")

ss <- stabsel(mod, q = 50)
}

### selected:
fm <- paste("bbs(CMEDV, boundary.knots = bk, df = 2.1) ~",
      " bbs(LSTAT, df = 2.1)" ,
      "+ bbs(CRIM, df = 2.1) ",
      "+ bbs(RM, df = 2.1)")
fm <- as.formula(fm)

mod <- ctm(fm, data = BH,
              family = Binomial(link = "probit"),  monotone = FALSE,
              control = boost_control(mstop = 250, trace = TRUE, nu = 0.2),
              ngrid = 25, weights = weights)

cv <- cvrisk(mod)

save(cv, file = "ex_BostonHousing.Rda")

library("lattice")
library("lattice")
trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))

pdf("ex_BostonHousing.pdf")
plot(mod, ylab = "Median Housing Value")
dev.off()
