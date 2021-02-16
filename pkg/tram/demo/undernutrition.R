### UNDERNUTRITION DATA DEMO
library("tram") 
library("mvtnorm")
library("colorspace")
library("latex2exp")

set.seed(42)
par(ask = TRUE)

### NOTE: simulated dataset. This is NOT the original
### undernutrition data. The original data can be obtained from
### https://dhsprogram.com/data/dataset/India_Standard-DHS_1999.cfm?flag=0
### for registered users.
### Preprocessing was performed as implemented in
###   system.file("india_preproc.R", package = "mboost")
### 
load(system.file("undernutrition.RData", package = "tram"))
summary(dat)
cageseq <- sort(unique(dat$cage))

################## MODELS ##################

### marginal models with linear shift in age
# m_stunting <- as.mlt(BoxCox(stunting2 ~ cage, data = dat, extrapolate = TRUE))
# m_wasting <- as.mlt(BoxCox(wasting2 ~ cage, data = dat, extrapolate = TRUE))
# m_underweight <- as.mlt(BoxCox(underweight2 ~ cage, data = dat, extrapolate = TRUE))

## marginal models with response-varying age effect
m_stunting <- as.mlt(BoxCox(stunting2 | cage ~ 1, data = dat, 
                            support = c(-4, 4), add = c(-2, 2)))
m_wasting <- as.mlt(BoxCox(wasting2 | cage ~ 1, data = dat, 
                           support = c(-4, 4), add = c(-2, 2)))
m_underweight <- as.mlt(BoxCox(underweight2 | cage ~ 1, data = dat, 
                               support = c(-4, 4), add = c(-2, 2)))

### parametrization for correlation coefficients
Bxlambda <- Bernstein_basis(numeric_var("cage", support = quantile(dat$cage, prob = c(.1, .9)),
                                        bounds = c(0, 100)), order = 6, extrapolate = TRUE)

### fitting joint model
m_full <- mmlt(m_stunting, m_wasting, m_underweight,
               formula = Bxlambda, data = dat, diag = FALSE)

### FAST ALTERNATIVE TO PARAMETRIC BOOTSTRAP
### sampling nsamp values from the asymptotic (normal) distribution of the parameters
nsamp <- 1000
V <- vcov(m_full)
V <- (V + t(V)) / 2
P <- rmvnorm(nsamp, mean = coef(m_full), sigma = V)
m_tmp <- m_full
CR <- vector(mode = "list", length = nrow(P))

nd <- data.frame(cage = cageseq)
ptm_npb <- system.time(
  for (i in 1:nsamp) {
    cf <- P[i, ]
    mi <- 1:length(m_full$pars$mpar)
    mcf <- cf[mi]
    vcf <- matrix(cf[-mi], nrow = nrow(m_full$pars$cpar))
    m_tmp$par <- cf
    m_tmp$pars <- list(mpar = mcf, cpar = vcf)
    CR[[i]] <- coef(m_tmp, newdata = nd, type = "Corr")
  }
)


### save estimated corr coef here
r12s <- r13s <- r23s <- matrix(NA, nrow = length(cageseq), ncol = nsamp) 
for(l in 1:nsamp) {
  r12s[, l] <- CR[[l]][, 1]
  r13s[, l] <- CR[[l]][, 2]
  r23s[, l] <- CR[[l]][, 3]
}
### save Spearman rhos here
rs12s <- 6*(asin(0.5*r12s))/pi
rs13s <- 6*(asin(0.5*r13s))/pi
rs23s <- 6*(asin(0.5*r23s))/pi
### estimated Spearman rhos
Cor_m_full <- coef(m_full, newdata = nd, type = "Corr")
rs12est <- 6*asin(0.5*Cor_m_full[, 1])/pi
rs13est <- 6*asin(0.5*Cor_m_full[, 2])/pi
rs23est <- 6*asin(0.5*Cor_m_full[, 3])/pi

################## PLOTS ##################
### only choose 1, 3, 6, 9, 12, 24 months
nd <- data.frame(cage = as.double(c(1, 3, 6, 9, 12, 24)))

### grids for distribution and density evaluation
q_stunting <- mkgrid(m_stunting, n = 100)[[1]]
q_wasting <- mkgrid(m_wasting, n = 100)[[1]]
q_underweight <- mkgrid(m_underweight, n = 100)[[1]]

### MARGINAL DISTRIBUTIONS
par(mfrow = c(1, 3), mar = c(5.5, 6.5, 3.5, 1.5) - 1)
d_stunting <- predict(m_full, newdata = nd, marginal = 1, 
                      type = "distribution", q = q_stunting)
col <- diverging_hcl(7, "Berlin")[-4]
plot(q_stunting, d_stunting[, 1], type = "n", ylim = c(0, 1), xlim = c(-5, 5),
     cex.lab = 2.5, cex.axis = 2,  panel.first = grid(),
     xlab = expression(paste(y[stunting], sep = "")),
     ylab = expression(paste("F(", y[stunting], "|age)", sep = "")))
for(i in 1:nrow(nd)) {
  lines(q_stunting, d_stunting[, i], col = col[i], lwd = 2)
}

d_wasting <- predict(m_full, newdata = nd, marginal = 2, 
                     type = "distribution", q = q_wasting)
plot(q_wasting, d_wasting[, 1], type = "n", ylim = c(0, 1), xlim = c(-5, 5),
     cex.lab = 2.5, cex.axis = 2,  panel.first = grid(),
     xlab = expression(paste(y[wasting], sep = "")),
     ylab = expression(paste("F(", y[wasting], "|age)", sep = "")))
for(i in 1:nrow(nd)) {
  lines(q_wasting, d_wasting[, i], col = col[i], lwd = 2)
}

d_underweight <- predict(m_full, newdata = nd, marginal = 3, 
                         type = "distribution", q = q_underweight)
plot(q_underweight, d_underweight[, 1], type = "n", ylim = c(0, 1), xlim = c(-5, 5),
     cex.lab = 2.5, cex.axis = 2,  panel.first = grid(),
     xlab = expression(paste(y[underweight], sep = "")),
     ylab = expression(paste("F(", y[underweight], "|age)", sep = "")))

for(i in 1:nrow(nd)) {
  lines(q_underweight, d_underweight[, i], col = col[i], lwd = 2)
}
legend("bottomright", legend = c(1, 3, 6, 9, 12, 24), title = "cage month", 
       col = col, bty = "n", lwd = 2, seg.len = .9, cex = 1.5)

### MARGINAL DENSITIES
par(mfrow = c(1, 3), mar = c(5.5, 6.5, 3.5, 1.5) - 1)
de_stunting <- predict(m_full, newdata = nd, marginal = 1, 
                       type = "density", q = q_stunting)
plot(q_stunting, de_stunting[, 1], type = "n", ylim = c(0, .4), xlim = c(-5, 5),
     cex.lab = 2.5, cex.axis = 2,  panel.first = grid(),
     xlab = expression(paste(y[stunting], sep = "")),
     ylab = expression(paste("f(", y[stunting], "|age)", sep = "")))
for(i in 1:nrow(nd)) {
  lines(q_stunting, de_stunting[, i], col = col[i], lwd = 2)
}

de_wasting <- predict(m_full, newdata = nd, marginal = 2, 
                      type = "density", q = q_wasting)
plot(q_wasting, de_wasting[, 1], type = "n", ylim = c(0, .4), xlim = c(-5, 5),
     cex.lab = 2.5, cex.axis = 2,  panel.first = grid(),
     xlab = expression(paste(y[wasting], sep = "")),
     ylab = expression(paste("f(", y[wasting], "|age)", sep = "")))
for(i in 1:nrow(nd)) {
  lines(q_wasting, de_wasting[, i], col = col[i], lwd = 2)
}

de_underweight <- predict(m_full, newdata = nd, marginal = 3, 
                          type = "density", q = q_underweight)
plot(q_underweight, de_underweight[, 1], type = "n", ylim = c(0, 0.4), xlim = c(-5, 5),
     cex.lab = 2.5, cex.axis = 2,  panel.first = grid(),
     xlab = expression(paste(y[underweight], sep = "")),
     ylab = expression(paste("f(", y[underweight], "|age)", sep = "")))
for(i in 1:nrow(nd)) {
  lines(q_underweight, de_underweight[, i], col = col[i], lwd = 2)
}
legend("topright", legend = c(1, 3, 6, 9, 12, 24), title = "cage month", 
       col = col, bty = "n", lwd = 2, seg.len = .9, cex = 1.5)

### correlation coefficients from NONPARAMETRIC BOOTSTRAP
par(mfrow = c(1, 3), mar = c(5.5, 7.9, 3.5, 1.5) - 1)
plot(cageseq, apply(rs12s, MARGIN = 1, FUN = "mean"), type = "l", lwd = 2, 
     xlab = "age",
     ylab = TeX('$\\rho_{stunting,wasting}^S(age)$'),
     cex.axis = 2.5, cex.lab = 2, cex.main = 2,
     ylim = c(-0.4, 0.15),
     panel.first = grid())
lines(cageseq, apply(rs12s, MARGIN = 1, FUN = function(x){quantile(x, prob = 0.025)}),
      lty = 2, lwd = 2)	
lines(cageseq, apply(rs12s, MARGIN = 1, FUN = function(x){quantile(x, prob = 0.975)}),
      lty = 2, lwd = 2)
# estimated Spearman coef from m_full
lines(cageseq, rs12est, lty = 2, lwd = 2, col = "red")

plot(cageseq, apply(rs13s, MARGIN = 1, FUN = "mean"), type = "l", lwd = 2, 
     xlab = "age",
     ylab = TeX('$\\rho_{stunting,underweight}^S(age)$'),
     cex.axis = 2.5, cex.lab = 2, cex.main = 2,
     ylim = c(0.45, 1.0),
     panel.first = grid())	
lines(cageseq, apply(rs13s, MARGIN = 1, FUN = function(x){quantile(x, prob = 0.025)}),
      lty = 2, lwd = 2)	
lines(cageseq, apply(rs13s, MARGIN = 1, FUN = function(x){quantile(x, prob = 0.975)}),
      lty = 2, lwd = 2)	
lines(cageseq, rs13est, lty = 2, lwd = 2, col = "red")

plot(cageseq, apply(rs23s, MARGIN = 1, FUN = "mean"), type = "l", lwd = 2, 
     xlab = "age",
     ylab = TeX('$\\rho_{wasting,underweight}^S(age)$'),
     cex.axis = 2.5, cex.lab = 2, cex.main = 2,
     ylim = c(0.4, 0.95),
     panel.first = grid())
lines(cageseq, apply(rs23s, MARGIN = 1, FUN = function(x){quantile(x, prob = 0.025)}),
      lty = 2, lwd = 2)	
lines(cageseq, apply(rs23s, MARGIN = 1, FUN = function(x){quantile(x, prob = 0.975)}),
      lty = 2, lwd = 2)	
lines(cageseq, rs23est, lty = 2, lwd = 2, col = "red")
legend("topleft", legend = c("95% CI", "mean", "ML estimate"), 
       col = c("black", "black", "red"), lty = c(2, 1, 2),
       bty = "n", lwd = c(2, 2, 2), cex = 1.5)

### warnings can be safely ignored

