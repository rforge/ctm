### 5D DEMO
library("mvtnorm")
library("VGAM")
library("parallel")
library("tram")
library("lattice")
library("latticeExtra")

################## DATA ##################
set.seed(29081975)

n <- 1000
J <- 5
repl <- 100

f21 <- function(x) x^2
f31 <- function(x) -x
f32 <- function(x) x^3 - x
f41 <- function(x) 0*x
f42 <- function(x) 0*x
f43 <- function(x) 0*x
f51 <- function(x) 0*x
f52 <- function(x) 0*x
f53 <- function(x) 0*x
f54 <- function(x) 0*x

### given x (lambda in paper), computes the interesting matrices
mats <- function(x){
  L <- matrix(c(1, f21(x), f31(x), f41(x), f51(x),
                0, 1, f32(x), f42(x), f52(x),
                0, 0, 1, f43(x), f53(x),
                0, 0, 0, 1, f54(x),
                0, 0, 0, 0, 1), nrow = J)
  Linv <- forwardsolve(L, diag(J))
  Sigma <- tcrossprod(Linv)
  Corr <- cov2cor(Sigma)
  return(list(L = L, Linv = Linv, Sigma = Sigma, Corr = Corr))
}

xseq <- seq(-0.9, 0.9, by = 0.1)
x <- sample(xseq, size = n, replace = TRUE)
mats_x <- lapply(x, mats)
a1 <- exp(2)
a2 <- exp(1.8)
a3 <- exp(1.5)
b1 <- exp(1)
b2 <- exp(0)
b3 <- exp(-0.9)
p1 <- exp(1.3)
p2 <- exp(0.9)
p3 <- exp(1)

data <- list()
for(i in 1:repl) {
  z <- tildez <- rmvnorm(n, mean = rep(0, J), sigma = diag(J))
  y1 <- y2 <- y3 <- y4 <- y5 <- rep(0, n)
  u1 <- u2 <- u3 <- u4 <- u5 <- rep(0, n)
  for(l in 1:n) {
    mats_l <- mats_x[[l]]
    tildez[l, ] <- mats_l$Linv %*% z[l, ]
    u1[l] <- pnorm(tildez[l, 1]) # ecdf(y1)(y1) - 1/(2*n)
    u2[l] <- pnorm(tildez[l, 2], mean = 0, sd = sqrt(mats_l$Sigma[2, 2])) # ecdf(y2)(y2) - 1/(2*n)
    u3[l] <- pnorm(tildez[l, 3], mean = 0, sd = sqrt(mats_l$Sigma[3, 3])) # ecdf(y3)(y3) - 1/(2*n)
    u4[l] <- pnorm(tildez[l, 4], mean = 0, sd = sqrt(mats_l$Sigma[4, 4]))
    u5[l] <- pnorm(tildez[l, 5], mean = 0, sd = sqrt(mats_l$Sigma[5, 5]))
    y1[l] <- qdagum(u1[l], shape1.a = a1, scale = b1, shape2.p = p1)
    y2[l] <- qdagum(u2[l], shape1.a = a2, scale = b2, shape2.p = p2)
    y3[l] <- qdagum(u3[l], shape1.a = a3, scale = b3, shape2.p = p3)
    y4[l] <- qgamma(u4[l], shape = 10.5, scale = 0.5)
    y5[l] <- qt(u5[l], 5)
  }
  data[[i]] <- data.frame(u1 = u1, u2 = u2, u3 = u3, u4 = u4, u5 = u5, 
                          y1 = y1, y2 = y2, y3 = y3, y4 = y4, y5 = y5, x = x)
}

################## MODELS ##################
run_mmlt <- function(i) {
  data_i <- data[[i]]
  
  ### Bernstein bases
  By <- lapply(c("y1","y2", "y3", "y4"), function(y) {
    v <- numeric_var(y, support = quantile(data_i[[y]], prob = c(.1, .9)),
                     bounds = c(0, Inf))
    Bernstein_basis(var = v, order = 6, ui = "increasing")
  })
  v <- numeric_var("y5", support = quantile(data_i[["y5"]], prob = c(.1, .9)),
                   bounds = c(-Inf, Inf))
  By5 <- Bernstein_basis(var = v, order = 6, ui = "increasing")
  By <- c(By, By5)
  Bx_shift <- Bernstein_basis(numeric_var("x", support = c(-.8, .8)), order = 6, 
                              ui = "zero")
  Bx_lambda <- Bernstein_basis(numeric_var("x", support = c(-.8, .8)), order = 6)
  
  ### marginal models
  ctm_y1 <- ctm(By[[1]], shift = Bx_shift, todistr = "Normal")
  m_y1 <- mlt(ctm_y1, data = data_i) 
  ctm_y2 <- ctm(By[[2]], shift = Bx_shift, todistr = "Normal")
  m_y2 <- mlt(ctm_y2, data = data_i)
  ctm_y3 <- ctm(By[[3]], shift = Bx_shift, todistr = "Normal")
  m_y3 <- mlt(ctm_y3, data = data_i)
  ctm_y4 <- ctm(By[[4]], shift = Bx_shift, todistr = "Normal")
  m_y4 <- mlt(ctm_y4, data = data_i)
  ctm_y5 <- ctm(By[[5]], shift = Bx_shift, todistr = "Normal")
  m_y5 <- mlt(ctm_y5, data = data_i)
  
  ### full model
  ptm <- system.time(
    m <- mmlt(m_y1, m_y2, m_y3, m_y4, m_y5, formula = Bx_lambda, data = data_i)
  )
  
  ### predict correlations
  rho <- coef(m, newdata = data.frame(x = xseq), type = "Lambda")
  ret <- list(rho = rho, ptm = ptm[3])
  return(ret)
}

resMLT <- mclapply(1:repl, FUN = run_mmlt, mc.cores = 4)



################## PLOTS ##################
Jp <- J*(J-1)/2

res_all <- data.frame(x = rep(xseq, repl + 1), 
                      repl = rep(1:(repl + 1), each = length(xseq)))
dd <- c()
for (i in 1:repl) {
  dd <- rbind(dd, resMLT[[i]]$rho)
}
dd <- rbind(dd, cbind(f21(xseq), f31(xseq), f32(xseq), 
                      matrix(0, nrow = length(xseq), ncol = (Jp - 3))))

id <- c()
for (i in 2:J) {
  for (k in 1:(i-1)) {
    id <- c(id, paste("l", i, k, sep = ""))
  }
}

res_all <- as.data.frame(c(res_all, as.data.frame(dd)))
colnames(res_all) <- c("x", "repl", id)

fm <- as.formula(paste(paste(id, collapse = "+"), "~ x " ))

# pdf("sim5d_lattice.pdf", height = 6, width = 12)
par(mar = c(5.5, 6.0, 3.5, 1.2) - 1)
panel_f <- function(x, y, repl = 100) {
  xseq <- seq(-0.9, 0.9, by = 0.1)
  panel.grid(h = -1, v = -1)
  for (i in 1:(repl+1)) {
    idx <- ((1+length(xseq)*(i-1)):(length(xseq)*i))
    panel.lines(x[idx], y[idx], col = "black", alpha = .2)
    if(i == (repl+1))
      panel.lines(x[idx], y[idx], col = "black", lwd = 2)
  }
}

xyplot(fm, group = repl, data = res_all, as.table = TRUE,
       between = list(x = .1), panel = panel_f, layout = c(5, 2),
       xlab = "x", ylab = expression(paste(lambda[ij], "(x)", sep = "")),
       strip = strip.custom(bg = "transparent",
                            factor.levels = 
                              c(expression(paste(lambda[21], "(x)", sep = "")),
                                expression(paste(lambda[31], "(x)", sep = "")),
                                expression(paste(lambda[32], "(x)", sep = "")),
                                expression(paste(lambda[41], "(x)", sep = "")),
                                expression(paste(lambda[42], "(x)", sep = "")),
                                expression(paste(lambda[43], "(x)", sep = "")),
                                expression(paste(lambda[51], "(x)", sep = "")),
                                expression(paste(lambda[52], "(x)", sep = "")),
                                expression(paste(lambda[53], "(x)", sep = "")),
                                expression(paste(lambda[54], "(x)", sep = "")))),
       scales = list(y = list(relation = "free"), x = list(relation = "free"),
                     alternating = 1))
# dev.off()

### boxplot 
logMSE <- vector(mode = "list", length = Jp)
for (i in 1:length(logMSE)) {
  logMSE[[i]] <- vector(repl, mode = "numeric")
}

logMSE[[1]] <- unlist(lapply(1:repl, FUN = function(x){
  log(sum((resMLT[[x]]$rho[, 1] - f21(xseq))^2/length(xseq)))}))
logMSE[[2]] <- unlist(lapply(1:repl, FUN = function(x){
  log(sum((resMLT[[x]]$rho[, 2] - f31(xseq))^2/length(xseq)))}))
logMSE[[3]] <- unlist(lapply(1:repl, FUN = function(x){
  log(sum((resMLT[[x]]$rho[, 3] - f32(xseq))^2/length(xseq)))}))
for(i in 4:length(logMSE)) {
  logMSE[[i]] <- unlist(lapply(1:repl, FUN = function(x){
    log(sum((resMLT[[x]]$rho[, i] - 0*(xseq))^2/length(xseq)))}))
}

### boxplot normal
# postscript("sim5dMSE.pdf", paper = "special", height = 5, width = 10)
par(mfrow = c(1, 1), mar = c(5.5, 6.0, 3.5, 1.2) - 1)
boxplot(sqrt(exp(as.matrix(as.data.frame(logMSE)))), ylim = c(0, 0.2), main = "RMSE",
        names = c(expression(paste(lambda[21], "(x)", sep = "")),
                  expression(paste(lambda[31], "(x)", sep = "")),
                  expression(paste(lambda[32], "(x)", sep = "")),
                  expression(paste(lambda[41], "(x)", sep = "")),
                  expression(paste(lambda[42], "(x)", sep = "")),
                  expression(paste(lambda[43], "(x)", sep = "")),
                  expression(paste(lambda[51], "(x)", sep = "")),
                  expression(paste(lambda[52], "(x)", sep = "")),
                  expression(paste(lambda[53], "(x)", sep = "")),
                  expression(paste(lambda[54], "(x)", sep = ""))),
        ylab = expression(paste("RMSE(", lambda[ij],"(x))", sep = "")),
        cex.axis = 1.75, cex.lab = 2, cex.main = 2)
# dev.off()