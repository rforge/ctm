library("cotram")
library("MASS")

## set-up plots
library("reshape2")
library("lattice")
library("colorspace")
col <- qualitative_hcl(3, c(240, 0), l = 60)

## settings
set.seed(29)
Nsim <- 100

## coefficients
b0 <- 1.2
b1 <- 0.8
theta <- 3

## transformation function dgp
pY <- ppois(0:100, lambda = exp(b0))
h.lo <- mlt:::.Logistic()[["q"]](pY) # logit link
h.cll <- mlt:::.MinExtrVal()[["q"]](pY) # cloglog link
h.ll <- mlt:::.MaxExtrVal()[["q"]](pY) # loglog link
h.pr <- mlt:::.Normal()[["q"]](pY) # probit link

## data-generating processes
dgp <- function(x = runif(1000, min = 0, max = 1)){
  log.mu <- b0 + b1 * x
  n <- length(x)
  
  ### Poisson - DGP
  y.p <- rpois(n = n, lambda = exp(log.mu))
  
  ### negbin - DGP
  y.nb <- rnbinom(n = n, mu = exp(log.mu), size = theta)
  
  ### logit - DGP
  h.lo_m <- matrix(h.lo, nrow = length(h.lo), ncol = length(x))
  p.lo <-  mlt:::.Logistic()[["p"]](t(h.lo_m) - b1 * x) # inv. logit link
  y.lo <- max.col(-(p.lo - runif(n))^2) - 1L
  
  mlo <- cotram(y.lo ~ x, data = data.frame(y.lo = y.lo, x = x), 
                method = "logit", fixed = c("x" = b1))
  
  ### cloglog - DGP
  h.cll_m <- matrix(h.cll, nrow = length(h.cll), ncol = length(x)) 
  p.cll <-  mlt:::.MinExtrVal()[["p"]](t(h.cll_m) - b1 * x) # inv. cloglog link
  y.cll <-  max.col(-(p.cll - runif(n))^2) - 1L
  
  mcll <- cotram(y.cll ~ x, data = data.frame(y.cll = y.cll, x = x), 
                 method = "cloglog", fixed = c("x" = b1))
  
  ### loglog - DGP
  h.ll_m <- matrix(h.ll, nrow = length(h.ll), ncol = length(x)) 
  p.ll <-  mlt:::.MaxExtrVal()[["p"]](t(h.ll_m) - b1 * x) # inv. loglog link
  y.ll <-  max.col(-(p.ll - runif(n))^2) - 1L
  
  mll <- cotram(y.ll ~ x, data = data.frame(y.ll = y.ll, x = x), 
                method = "loglog", fixed = c("x" = b1))
  
  ### probit - DGP
  h.pr_m <- matrix(h.pr, nrow = length(h.pr), ncol = length(x)) 
  p.pr <-  mlt:::.Normal()[["p"]](t(h.pr_m) - b1 * x) # inv. probit link
  y.pr <- max.col(-(p.pr - runif(n))^2) - 1L
  
  mpr <- cotram(y.pr ~ x, data = data.frame(y.pr = y.pr, x = x), 
                method = "probit", fixed = c("x" = b1))
  
  ### data frame
  ret <- data.frame(x = x, mu = exp(log.mu),
                    y.p = y.p, y.nb = y.nb, y.lo = y.lo,
                    y.cll = y.cll, y.ll = y.ll, y.pr = y.pr)
  
  attr(ret, "mlo") <- mlo
  attr(ret, "mcll") <- mcll
  attr(ret, "mll") <- mll
  attr(ret, "mpr") <- mpr
  return(ret)
}

## model out-of-sample log-likelihood
logLikFUN <- function(model, lhs, newdata){ 
  
  if("cotram" %in% class(model)){
    val <- logLik(model, newdata = newdata)
    
  } else if(grepl("Negative Binomial", family(model)$family)){
    val <- sum(dnbinom(newdata[,lhs],
                       mu = predict(model, newdata = newdata, type = "response"),
                       size = model$theta, log = TRUE))
    
  } else if(family(model)$family %in% "poisson"){
    val <- sum(dpois(newdata[,lhs],
                     lambda = predict(model, newdata = newdata, type = "response"),
                     log = TRUE))
  }
  return(val)}

## out-of-sample log-likelihood of generating process
logLik_refFUN <- function(dgp, newdata){
  if(dgp == "lo"){
    mlo <- attr(newdata, "mlo")
    val <- logLik(mlo, newdata = newdata, parm = coef(as.mlt(mlo), fixed = TRUE))
    
  }else if(dgp == "cll"){
    mcll <- attr(newdata, "mcll")
    val <- logLik(mcll, newdata = newdata, parm = coef(as.mlt(mcll), fixed = TRUE))
    
  }else if(dgp == "ll"){
    mll <- attr(newdata, "mll")
    val <- logLik(mll, newdata = newdata, parm = coef(as.mlt(mll), fixed = TRUE))
    
  }else if(dgp == "pr"){
    mpr <- attr(newdata, "mpr")
    val <- logLik(mpr, newdata = newdata, parm = coef(as.mlt(mpr), fixed = TRUE))
    
  }else if(dgp == "nb"){
    val <- sum(dnbinom(newdata$y.nb, mu = newdata$mu, size = theta, log = TRUE))
    
  }else if(dgp == "p"){
    val <- sum(dpois(newdata$y.p, lambda = newdata$mu, log = TRUE))
  }
  return(val)
}

## set-up
dgps <- c("p", "nb", "lo", "cll", "ll", "pr")
mods <- paste0("m", dgps)

logLiks <- vector(mode = "list", length = length(dgps))
names(logLiks) <- dgps
for (d in dgps){
  logLiks[[d]] <- as.data.frame(matrix(NA, nrow = Nsim, ncol = length(dgps),
                                   dimnames = list(1:Nsim, mods)))
  logLiks[[d]]$dgp <- NA
}

## run Nsim-times
for (i in c(1:Nsim)){
  print(i)
  
  df <- dgp()
  
  ### split into training and validation data-set
  trainID <- sample(1:nrow(df), size = 0.25 * nrow(df))
  df_train <- df[trainID,]
  df_test <- df[-trainID,]
  
  ### for each data-generating process
  for (d in dgps){
    
    ### model formula
    lhs <- paste0("y.", d)
    fm <-  as.formula(paste(lhs, "~ x"))
    
    ### count regression models
    mp <- glm(fm, data = df_train, family = poisson(link = "log"))
    mnb <- glm.nb(fm, data = df_train, link = "log")
    
    mlo <- cotram(fm, data = df_train, method = "logit")
    mcll <- cotram(fm, data = df_train, method = "cloglog")
    mll <- cotram(fm, data = df_train, method = "loglog")
    mpr <- cotram(fm, data = df_train, method = "probit")
    
    fit <- list(mp = mp, mnb = mnb, mlo = mlo, mcll = mcll, mll = mll, mpr = mpr)
    
    ### out-of-sample log-likelihood of models
    logLiks.d <- sapply(fit, logLikFUN, lhs = lhs, newdata = df_test)
    
    ### centered out-of-sample log-likelihood
    logLiks_diff.d <- logLiks.d - logLik_refFUN(newdata = df_test, dgp = d)
    
    logLiks[[d]][i, names(logLiks_diff.d)] <- logLiks_diff.d
    logLiks[[d]][i, "dgp"] <- d
  }
}


## plots of the out-of-sample log-likelihoods
tmp_logLiks <- do.call("rbind", logLiks)
tmp_logLiks$id <- rownames(tmp_logLiks)
tmp_logLiks <- melt(tmp_logLiks, id.vars = c("id", "dgp"), value.name = "logLik", variable.name = "fit")

tmp_logLiks$dgp <- ordered(as.character(tmp_logLiks$dgp),
                   levels = c("ll", "pr", "lo", "p", "nb", "cll"),
                   labels = paste(c("loglog", "probit", "logit", "Poisson", "neg. binomial", "cloglog"), 
                                  "DGP", sep = " - ")
)
tmp_logLiks$fit <- ordered(as.character(tmp_logLiks$fit),
                   levels = c("mp", "mnb", "mlo", "mcll", "mll", "mpr"),
                   labels = c(paste(c("Poisson", "neg. binomial"), "model"),
                              paste(c("logit", "cloglog", "loglog", "probit"),
                                    "transformation model")))
bwplot(logLik  ~ fit | dgp, data = tmp_logLiks, groups = id, type = "l", lwd = 2,
       ylab = "Centered out-of-sample log-likelihood",
       panel = function(x, y, groups, subscripts, ...) {
         panel.abline(h = 0, col = "gray70", lty = 3)
         panel.bwplot(x = x, y = y, ...)
         tapply(1:length(y), groups[subscripts], function(i) {
           llines(
             x = 1:nlevels(x),
             y = y[i][order(x[i])],
             col = rgb(.1, .1, .1, .1))})},
       layout = c(3, 2),
       scales = list(x = list(rot = 45), y = list(relation = 'free')))



## conditional distribution of the dgps ##
x <-  seq(0, 1, length.out = 50)
q <- 0:20

df.dgps <- vector(mode = "list", length = length(dgps))
names(df.dgps) <- dgps

for (d in dgps){
  nd <- expand.grid(q = q, x = x)
  nd$log.mu <- b0 + b1 * nd$x
  
  if (d == "p") nd$p <- ppois(nd[, "q"], lambda = exp(nd$log.mu))
  if (d == "nb") nd$p <- pnbinom(nd[, "q"], mu = exp(nd$log.mu), size = theta)
  if (d %in% c("lo", "cll", "ll", "pr")) {
    m  <- attr(df, paste0("m", d))
    nd$p <- mapply(FUN <- function(x, q){
      c(predict(m, newdata = data.frame(x = x), 
                type = "distribution", q = q))},
      x = nd$x, q = nd$q)
  }
  nd$dgp <- d
  df.dgps[[d]] <- as.data.frame(nd)
}


## plots of cdf of dgps
tmp_dgps <- do.call(rbind, df.dgps)
tmp_dgps$dgp <- ordered(as.character(tmp_dgps$dgp),
                   levels = c("ll", "pr", "lo", "p", "nb", "cll"),
                   labels = paste(c("loglog", "probit", "logit", "Poisson", "neg. binomial", "cloglog"), 
                                  "DGP", sep = " - "))

wireframe(p ~ q * x | dgp, data = tmp_dgps,
          xlab = list(label = "Counts", cex = .75, rot = 30),
          ylab = list(label = "x", cex = .75, rot = -30),
          drape = TRUE,
          col.regions = sequential_hcl(200, "Teal"),
          col = "transparent",
          scales = list(arrows = FALSE,
                        cex = .5,
                        x = list(tck = .8),
                        y = list(tck = .8, at = seq(0, 1, .2), label = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")),
                        z = list(tck = .8, distance = 5)
          ), layout = c(3, 2), 
          ylab.right = "Probability",
          par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                   ylab.right = 3)))


warnings()
