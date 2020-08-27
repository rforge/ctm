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


##### Estimation with BayesX
# ### pathwd has to be adapted depending on the user
# pathwd <- "/home/luisa/multivariate/mlt_multivariate/RevisionSJS/codes/simulation/2d"
# ## path where bayesx is stored (you need a linux self-compiled SVN version of BayesX)
# pathbayesx <- "/home/luisa/bayesx/bayesx/trunk/bayesx"
# ## path for prg's to be run in bayesx
# pathprg <- paste(pathwd, "/prg/", sep = "")
# ## path where bayesx results should be stored
# pathresults <- paste(pathwd, "/resultsBayesX/", sep = "")
# 
# source(paste(pathwd, "utilsBayesX.r", sep = "/"))
# 
# library("BayesX")
# family <- "gaussian"
# fm_b <- "const"
# fm_b2 <- "const + x(pspline, lambda = 100)"
# 
# runBayesX <- function(i)  {
#   outdata <- list()
#   batchname <- datname <- paste("rep", i, sep = "")
#   dattmp <- data[[i]]
#   dattmp$w <- 1
#   datpred <- data.frame(u1 = 0, u2 = 0, y1 = 10, y2 = 10, x = xseq, w = 0)
#   dattmp <- rbind(dattmp, datpred)
#   write.table(dattmp, paste(pathprg, batchname, ".raw", sep = ""), quote = FALSE, row.names = FALSE)
#   filesBayesX(i, pathprg, datname, fm_b, fm_b2, batchname, family)
#   
#   ptm <- system.time(
#     system(paste(pathbayesx, paste(pathprg, batchname, ".txt", sep = "")))
#   )
#   outdata$ptm <- ptm
#   
#   if(file.exists(paste(pathresults, batchname, "_MAIN_rho_REGRESSION_y", "1_predict.res", sep = ""))) {
#     pred <- read.table(paste(pathresults, batchname, "_MAIN_rho_REGRESSION_y", "1_predict.res", sep = ""),
#                        header = TRUE)
#     pred <- pred[pred$w == 0, ]
#     outdata$pred <- pred
#     outdata$rho <- pred$pmean_param_rho
#     outdata$a <- -outdata$rho/sqrt(1 - outdata$rho^2)
#     outdata$x <- pred$x
#     pred <- read.table(paste(pathresults, batchname, "_MAIN_rho_REGRESSION_y",
#                              "1_nonlinear_pspline_effect_of_x_param.res", sep = ""),
#                        header = TRUE)
#     source(paste(pathresults, batchname,
#                  "_MAIN_rho_REGRESSION_y1_nonlinear_pspline_effect_of_x_basisR.res",
#                  sep = ""))
#     Bseq <- BayesX.design.matrix(xseq)
#     frho <- Bseq %*% pred$pmean
#     outdata$frho <- frho - mean(frho)
#     pred <- read.table(paste(pathresults, batchname, "_MAIN_rho_REGRESSION_y",
#                              "1_LinearEffects.res", sep = ""), header = TRUE)
#     outdata$constrho <- pred
#     outdata$predictor <- frho + pred$pmean
#   } else {
#     outdata$pred <- NA
#   }
#   save(outdata, file = paste(pathresults, "res_rep_", i, ".RData", sep = ""))
#   return(outdata)
# }
# 
# resBayesX <- mclapply(1:repl, FUN = runBayesX, mc.cores = 4)
# # save(resBayesX, file = "resBayesX.RData")
# 


##### Estimation with VGAM
# set.seed(29081975)
# runVGAM <- function(i) {
#   dattmp <- data[[i]]
#   ptm <- proc.time()
#   ret1 <- vglm(formula = y1 ~ 1, family = dagum, data = dattmp, half.stepsizing = TRUE,
#                stepsize = 0.5, maxit = 300, coefstart = c(0, 0, 0), noWarning = TRUE)
#   param1 <- exp(coef(ret1))
#   F1 <- pdagum(dattmp$y1, scale = param1[1], shape1.a = param1[2], shape2.p = param1[3],
#                lower.tail = TRUE, log.p = FALSE)
#   dattmp$F1 <- F1
#   clist <- list("(Intercept)" = diag(3), "sm.ps(x)" = rbind(1, 0, 0))
#   ret2 <- vglm(y2 ~ sm.ps(x), family = dagum, data = dattmp, constraints = clist,
#                epsilon = 1e-04, maxit = 300, checkwz = FALSE, half.stepsizing = TRUE,  
#                stepsize = stepsizes[i], trace = TRUE, Maxit.outer = 200)
#   
#   param2 <- exp(coef(ret2)[c(2,3)])
#   preddata <- data.frame(y2 = 0, x = xseq)
#   pred <- data.frame(x = xseq, predictor = predict(ret2, type = "link", newdata = preddata))
#   pred$b2 <- exp(pred$predictor.loglink.scale.)
#   F2 <- rep(0, length(dattmp$x))
#   for(kk in 1:length(dattmp$x)) {
#     findb <- pred$b2[which(xseq == dattmp$x[kk])]
#     F2[kk] <- pdagum(dattmp$y2[kk], scale = findb, shape1.a = param2[1],
#                      shape2.p = param2[2], lower.tail = TRUE, log.p = FALSE)
#   }
#   dattmp$F2 <- F2
#   res <- vgam(formula = cbind(F1, F2) ~ sm.ps(x, ps.int = 17), family = binormalcop,
#               data = dattmp, checkwz = F, stepsize = 0.1, maxit = 100, Maxit.outer = 100,
#               noWarning = TRUE)
#   ptm <- proc.time() - ptm
#   
#   preddata <- data.frame(F1 = 0, F2 = 0, x = xseq)
#   pred <- data.frame(x = xseq, predictor = predict(res, type = "link", newdata = preddata))
#   pred$rho <- rhobitlink(pred$predictor, inverse = TRUE)
#   pred$a <- -pred$rho / sqrt(1 - pred$rho^2)
#   pred$ptm <- ptm[3]
#   return(pred)
# }
# # defining special stepsizes different than 0.15
# stepsizes <- rep(0.15, 100)
# stepsizes[44] <- 0.17
# stepsizes[78] <- 0.13
# stepsizes[89] <- 0.17
# stepsizes[98] <- 0.19
# 
# ### there are 2 replications for which we couldn't tune VGAM properly. we exclude them
# resVGAM <- mclapply(1:79, FUN = runVGAM, mc.cores = 4)
# resVGAM <- c(resVGAM, mclapply(81:90, FUN = runVGAM, mc.cores = 4))
# resVGAM <- c(resVGAM, mclapply(92:100, FUN = runVGAM, mc.cores = 4))
# 
# ### making up for 2 missing fits so that there is no influence on median times
# resVGAM <- c(resVGAM, mclapply(c(1, 98), FUN = runVGAM, mc.cores = 4))
# resVGAM[[99]]$ptm <- 0
# resVGAM[[100]]$ptm <- 10
# # save(resVGAM, file = "resVGAM.RData")



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