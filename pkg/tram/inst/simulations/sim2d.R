### 2D DEMO
library("mvtnorm")
library("VGAM")
library("parallel")
library("tram")

################## DATA ##################
set.seed(29081975)

n <- 1000
J <- 2
repl <- 100

### given lambda coefficient in Lambda, computes the interesting matrices
mats <- function(lambda){
  L <- matrix(c(1, lambda, 0, 1), nrow = J)
  Linv <- matrix(c(1, -lambda, 0, 1), nrow = J)
  Sigma <- matrix(c(1, -lambda, -lambda, 1+lambda^2), nrow = J)
  Corr <- matrix(c(1, -lambda/sqrt(1+lambda^2), -lambda/sqrt(1+lambda^2), 1), nrow = J)
  return(list(L = L, Linv = Linv, Sigma = Sigma, Corr = Corr))
}

xseq <- seq(-0.9, 0.9, by = 0.1)
x <- sample(xseq, size = n, replace = TRUE)
mats_x <- lapply(x^2, mats)
a1 <- exp(2)
a2 <- exp(1.8)
b1 <- exp(1)
b2 <- exp(0)
p1 <- exp(1.3)
p2 <- exp(0.9)

data <- list()
for(i in 1:repl) {
  z <- tildez <- rmvnorm(n, mean = rep(0, J), sigma = diag(J))
  y1 <- y2 <- rep(0, n)
  u1 <- u2 <- rep(0, n)
  for(l in 1:n) {
    mats_l <- mats_x[[l]]
    tildez[l, ] <- mats_l$Linv %*% z[l, ]
    u1[l] <- pnorm(tildez[l, 1]) # ecdf(y1)(y1) - 1/(2*n)
    u2[l] <- pnorm(tildez[l, 2], mean = 0, sd = sqrt(mats_l$Sigma[2, 2])) # ecdf(y2)(y2) - 1/(2*n)
    y1[l] <- qdagum(u1[l], shape1.a = a1, scale = b1, shape2.p = p1)
    y2[l] <- qdagum(u2[l], shape1.a = a2, scale = b2, shape2.p = p2)
  }
  data[[i]] <- data.frame(u1 = u1, u2 = u2, y1 = y1, y2 = y2, x = x)
}

################## MODELS ##################
run_mmlt <- function(i) {
  data_i <- data[[i]]
  
  ### Bernstein bases
  By <- lapply(c("y1","y2"), function(y) {
    v <- numeric_var(y, support = quantile(data_i[[y]], prob = c(.1, .9)),
                     bounds = c(0, Inf))
    Bernstein_basis(var = v, order = 6, ui = "increasing")  # change order?
  })
  Bx_shift <- Bernstein_basis(numeric_var("x", support = c(-.8, .8)), order = 6, ui = "zero")
  Bx_lambda <- Bernstein_basis(numeric_var("x", support = c(-.8, .8)), order = 6)
  
  ### marginal models
  ctm_y1 <- ctm(By[[1]], shift = Bx_shift, todistr = "Normal")
  m_y1 <- mlt(ctm_y1, data = data_i) 
  ctm_y2 <- ctm(By[[2]], shift = Bx_shift, todistr = "Normal")
  m_y2 <- mlt(ctm_y2, data = data_i)
  
  ### full model
  ptm <- system.time(m <- mmlt(m_y1, m_y2, formula = Bx_lambda, data = data_i))
  
  ### predict correlations
  rho <- coef(m, newdata = data.frame(x = xseq), type = "Lambda")
  ret <- list(rho = rho, ptm = ptm)
  return(ret)
}

system.time(
  resMLT <- mclapply(1:repl, FUN = run_mmlt, mc.cores = 4)
)


################## PLOTS ##################
postscript("sim2d.eps", paper = "special", height = 4, width = 12)
par(mfrow = c(1, 1), mar = c(5.5, 6.0, 3.5, 1.2) - 1)
plot(xseq, xseq^2, type = "n", ylim = c(-0.2, 1), main = "MCTM", 
     xlab = "x", ylab = expression(paste(lambda,"(x)", sep = "")),
     cex.axis = 1.75, cex.lab = 2, cex.main = 2)	
for(i in 1:repl) {
  lines(xseq, resMLT[[i]]$rho, col = rgb(.1, .1, .1, .1))
}
lines(xseq, xseq^2, lwd = 2)
dev.off()

logMSE <- list(MCTM = vector(repl, mode = "numeric"))
logMSE$MCTM <- unlist(lapply(1:repl, FUN = function(x){log(sum((resMLT[[x]]$rho - xseq^2)^2/
                                                                 length(xseq)))}))

postscript("sim2dMSE.eps", paper = "special", height = 6, width = 10)
par(mfrow = c(1, 1), mar = c(4.5, 6.0, 3.5, 1.2) - 1)
logMSE <- as.matrix(as.data.frame(logMSE))	
boxplot(sqrt(exp(logMSE)),
        main = "RMSE", ylab = expression(paste("RMSE(", lambda, "(x))", sep = "")),
        cex.axis = 1.75, cex.lab = 2, cex.main = 2)
dev.off()

postscript(paste(pathwd, "sim2dMSE2.eps", sep = ""), paper = "special", height = 10, width = 20)
par(mfrow = c(1, 2), mar = c(4.5, 6.0, 3.5, 1.2) - 1)
logMSE <- as.matrix(as.data.frame(logMSE))	
boxplot(exp(logMSE), main = "MSE", ylab = "MSE",
        cex.axis = 1.75, cex.lab = 2, cex.main=2)
boxplot((logMSE), main = "log(MSE)", ylab = "log(MSE)",
        cex.axis = 1.75, cex.lab = 2, cex.main=2)
dev.off()
