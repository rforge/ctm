## ----setopts, echo=FALSE, message=FALSE----------------------------------
knitr::opts_chunk$set(size = "footnotesize",
                      fig.width = 6, fig.height = 4, fig.align = "center")
knitr::opts_knit$set(global.par = TRUE)
options(digits = 6)

## ----defs, include=FALSE-------------------------------------------------
mycolors <- function(nr, type = "line") {
  cols <- list()
  cols[[1]] <- c(red = 0, green = 84, blue = 150)
  cols[[2]] <- c(red = 202, green = 108, blue = 24)
  out <- as.list(cols[[nr]])
  out$alpha <- switch(type, line = 255L, fill = 140L)
  out$maxColorValue <- 255
  do.call("rgb", out)
}

## ------------------------------------------------------------------------
data("sleepstudy", package = "lme4")
## plot

## ----message=FALSE-------------------------------------------------------
library("tramME")
sleep_lm <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
logLik(sleep_lm)

## ----echo=FALSE, fig.width=6, fig.height=3.5-----------------------------
par(mfrow = c(1, 2), cex = 0.75)
plot(sleep_lm, type = "trafo")
plot(sleep_lm, type = "distribution")

## ----message=FALSE-------------------------------------------------------
library("lme4")
sleep_lmer <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
logLik(sleep_lmer)

## ------------------------------------------------------------------------
cbind(coef = coef(sleep_lm, as.lm = TRUE),
      se = sqrt(diag(vcov(sleep_lm, as.lm = TRUE, pargroup = "fixef"))))

## ------------------------------------------------------------------------
summary(sleep_lmer)$coefficients

## ------------------------------------------------------------------------
VarCorr(sleep_lm, as.lm = TRUE) ## random effects
sigma(sleep_lm) ## residual SD
VarCorr(sleep_lmer)

## ------------------------------------------------------------------------
library("survival")
ub <- ceiling(sleepstudy$Reaction / 50) * 50
lb <- floor(sleepstudy$Reaction / 50) * 50
lb[ub == 200] <- 0
sleepstudy$Reaction_ic <- Surv(lb, ub, type = "interval2")
head(sleepstudy$Reaction_ic)

## ------------------------------------------------------------------------
sleep_lm2 <- LmME(Reaction_ic ~ Days + (Days | Subject), data = sleepstudy)
logLik(sleep_lm2)

## ------------------------------------------------------------------------
cbind(coef = coef(sleep_lm2, as.lm = TRUE),
      se = sqrt(diag(vcov(sleep_lm2, as.lm = TRUE, pargroup = "fixef"))))
sigma(sleep_lm2)
VarCorr(sleep_lm2, as.lm = TRUE)

## ------------------------------------------------------------------------
sleep_bc <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
logLik(sleep_bc)

## ----echo=FALSE,fig.width=4.5, fig.height=4------------------------------
par(mfrow = c(1, 1), cex = 0.8)
plot(sleep_bc, newdata = data.frame(Days = 0, Subject = 1),
     ranef = c(0, 0), type = "trafo", col = 1, lwd = 2,
     xlab = "Average reaction time (ms)", ylab = expression(h(y)))
plot(sleep_lm, newdata = data.frame(Days = 0, Subject = 1),
     ranef = c(0, 0), type = "trafo", add = TRUE,
     col = 1, lty = 2, lwd = 2)
rug(sleepstudy$Reaction, col = rgb(.1, .1, .1, .1))
legend("topleft", c("BoxCoxME", "LmME"), lty = c(1, 2), lwd = 2,
       bty = "n", cex = 0.8)

## ----eval=FALSE----------------------------------------------------------
#  ndraws <- 1000
#  nd <- expand.grid(
#    Reaction = seq(min(sleepstudy$Reaction), max(sleepstudy$Reaction), length.out = 100),
#    Days = 0:9,
#    Subject = 1)
#  
#  re <- simulate(sleep_bc, newdata = nd, nsim = ndraws, what = "ranef", seed = 100)
#  cp <- parallel::mclapply(re, function(x) {
#    predict(sleep_bc, newdata = nd, ranef = x, type = "distribution")
#  }, mc.cores = 8)
#  cp <- array(unlist(cp), dim = c(100, 10, ndraws))
#  mp_bc <- apply(cp, c(1, 2), mean)

## ----echo=FALSE, eval=FALSE----------------------------------------------
#  ## Numerical integration of the Box-Cox and Lm models. Does not run.
#  ndraws <- 1000
#  nd <- expand.grid(
#    Reaction = seq(min(sleepstudy$Reaction), max(sleepstudy$Reaction), length.out = 100),
#    Days = 0:9,
#    Subject = 1)
#  
#  ## --- Box-Cox
#  re <- simulate(sleep_bc, newdata = nd, nsim = ndraws, what = "ranef", seed = 100)
#  cp <- parallel::mclapply(re, function(x) {
#    predict(sleep_bc, newdata = nd, ranef = x, type = "distribution")
#  }, mc.cores = 8)
#  cp <- array(unlist(cp), dim = c(100, 10, ndraws))
#  mp_bc <- apply(cp, c(1, 2), mean)
#  
#  ## --- Lm
#  re <- simulate(sleep_lm, newdata = nd, nsim = ndraws, what = "ranef", seed = 100)
#  cp <- parallel::mclapply(re, function(x) {
#    predict(sleep_lm, newdata = nd, ranef = x, type = "distribution")
#  }, mc.cores = 8)
#  cp <- array(unlist(cp), dim = c(100, 10, ndraws))
#  mp_lm <- apply(cp, c(1, 2), mean)
#  
#  save(nd, mp_bc, mp_lm, file = "tramME/inst/vignette_data/marg_int.rda")

## ----echo=FALSE, message=FALSE, fig.width=7, fig.height=4.5--------------
load(system.file("vignette_data", "marg_int.rda", package = "tramME"))
dat <- nd
dat$Days <- paste0("Days = ", dat$Days)
dat$mp_lm <- c(mp_lm)
dat$mp_bc <- c(mp_bc)
dat2 <- sleepstudy
dat2$Days <- paste0("Days = ", dat2$Days)
library("ggplot2")
ggplot(dat, aes(x = Reaction)) +
  facet_wrap(~ Days) +
  geom_line(aes(y = mp_bc, colour = "BoxCoxME")) +
  geom_line(aes(y = mp_lm, colour = "LmME")) +
  stat_ecdf(aes(x = Reaction, colour = "ECDF"), data = dat2, geom = "step") +
  xlab("Average reaction time (ms)") +
  ylab("Marginal distribution") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.5, 0.3),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  scale_color_manual(
    values = c(rgb(0, 84, 150, maxColorValue = 255),
               rgb(.5, .5, .5, .5),
               rgb(202, 108, 24, maxColorValue = 255)),
    breaks = c("BoxCoxME", "LmME", "ECDF"))

