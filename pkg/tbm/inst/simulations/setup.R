
load("dgp.rda")
library("parallel")
RNGkind("L'Ecuyer-CMRG")

set.seed(290875)

tNOBS <- c(75, 150, 300)
NSIM <- 100
tPNON <- c(0, 5, 25)

tTD <- c("normal", "logistic", "minextrval")
tOR <- 1:2

MC_CORES <- 32
