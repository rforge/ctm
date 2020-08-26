### 2D DEMO
library("mvtnorm")
library("VGAM")
library("parallel")
library("tram")
library("lattice")
library("latticeExtra")

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
run_mmlt <- function(i, o_marginal, o_lambda) {
  data_i <- data[[i]]
  
  ### Bernstein bases
  By <- lapply(c("y1","y2"), function(y) {
    v <- numeric_var(y, support = quantile(data_i[[y]], prob = c(.1, .9)),
                     bounds = c(0, Inf))
    Bernstein_basis(var = v, order = o_marginal, ui = "increasing")
  })
  Bx_shift <- Bernstein_basis(numeric_var("x", support = c(-.8, .8)), order = 6, ui = "zero")
  Bx_lambda <- Bernstein_basis(numeric_var("x", support = c(-.8, .8)), order = o_lambda)
  
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

resMLT_6_6 <- mclapply(1:repl, FUN = run_mmlt, o_marginal = 6, o_lambda = 6, mc.cores = 4)
resMLT_6_3 <- mclapply(1:repl, FUN = run_mmlt, o_marginal = 6, o_lambda = 3, mc.cores = 4)
resMLT_12_3 <- mclapply(1:repl, FUN = run_mmlt, o_marginal = 12, o_lambda = 3, mc.cores = 4)


################## PLOTS ##################
tt1 <- tt2 <- tt3 <- c()
for(i in 1:repl) {
  tt1 <- c(tt1, resMLT_6_6[[i]]$rho)
  tt2 <- c(tt2, resMLT_6_3[[i]]$rho)
  tt3 <- c(tt3, resMLT_6_3[[i]]$rho)
}
res_all <- data.frame(x = rep(xseq, repl), repl = rep(1:100, each = length(xseq)))
res_all$MCTM66 <- tt1
res_all$MCTM63 <- tt2
res_all$MCTM123 <- tt3


# pdf("sim2d.pdf", paper = "special", height = 4, width = 12)
par(mar = c(5.5, 6.0, 3.5, 1.2) - 1)
panel_f <- function(x, y, repl = 100) {
  xseq <- seq(-0.9, 0.9, by = 0.1)
  panel.grid(h = -1, v = -1)
  for (i in 1:repl) {
    idx <- ((1+length(xseq)*(i-1)):(length(xseq)*i))
    panel.lines(x[idx], y[idx], col = "black", alpha = .2)
  }
  panel.lines(xseq, xseq^2, col = "black", lwd = 2)
}

xyplot(MCTM66 + MCTM63 + MCTM123 ~ x, group = repl,  data = res_all, outer = TRUE,
       between = list(x = 1), layout = c(3, 1), panel = panel_f,
       ylim = c(-0.2, 1), scales = list(x = list(relation = "free"),
                                        y = list(rot = 90)),
       xlab = "x", ylab = expression(paste(lambda, "(x)", sep = "")),
       strip = strip.custom(bg = "transparent",
                            factor.levels = c("MCTM-6/6", "MCTM-6/3", "MCTM-12/3")))
# dev.off()

logMSE <- list(MCTM66 = vector(repl, mode = "numeric"), 
               MCTM63 = vector(repl, mode = "numeric"), 
               MCTM123 = vector(repl, mode = "numeric"))

logMSE$MCTM66 <- unlist(lapply(1:repl, FUN = function(x){log(sum((resMLT_6_6[[x]]$rho - xseq^2)^2/
                                                                   length(xseq)))}))
logMSE$MCTM63 <- unlist(lapply(1:repl, FUN = function(x){log(sum((resMLT_6_3[[x]]$rho - xseq^2)^2/
                                                                   length(xseq)))}))
logMSE$MCTM123 <- unlist(lapply(1:repl, FUN = function(x){log(sum((resMLT_12_3[[x]]$rho - xseq^2)^2/
                                                                    length(xseq)))}))

# postscript("sim2dMSE_ext.pdf", paper = "special", height = 6, width = 11)
par(mfrow = c(1, 1), mar = c(4.5, 6.0, 3.5, 1.2) - 1)
logMSE <- as.matrix(as.data.frame(logMSE))	
boxplot(sqrt(exp(logMSE)), ylim = c(0, 0.2),
        main = "RMSE", ylab = expression(paste("RMSE(", lambda, "(x))", sep = "")),
        names = c("MCTM-6/6", "MCTM-6/3", "MCTM-12/3"), 
        cex.axis = 1.75, cex.lab = 2, cex.main = 2)
# dev.off()