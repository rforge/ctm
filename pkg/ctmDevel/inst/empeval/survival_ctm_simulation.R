
##********************************************************************
## CTM for survival data: Simulation 
##********************************************************************

## The R code given below enables the reader to reproduce the simulation results ## given in Section 3 of the main paper.
## All analysis were carried out in R version 3.0.2. Furthermore, one needs 
## to install the add-on packages ctmDevel (depending on mboostDevel), 
## survival (for the estimation of Cox models), prodlim (for the calculation of
## Kaplan-Meier estimators), and VGAM (for the density, quantile, and probabilty
## function of the Gumbel distribution).
 
###################################################
### code chunk number 1: Binomial_gumbel
###################################################

## The function Binomial_Gumbel defines a Family-object which is  compatible with
## the families from the mboost-package. The Binomial_gumbel-family is used for
## binary response variables and the corresponding link function is set to the
## distribution function of the Gumbel distribution:

library("VGAM")
Binomial_gumbel <- function(){
  biny <- function(y){
          if(!is.factor(y)) 
            stop("response is not a factor but ", sQuote("family = Binomial()"))
          if(nlevels(y) != 2) 
            stop("response is not a factor at two levels but ", 
                sQuote("family = Binomial()"))
    return(c(-1, 1)[as.integer(y)])
  }
  trf <- function(f) {
    pmax(pmin(f, 15), -2.5)
  }
  return(Family(
      ngradient = function(y, f, w = 1){
          trf <- function(f){
                 pmax(pmin(f, 15), -2.5)
          }
          y <- (y + 1)/2
          p <- pgumbel(trf(f))
          d <- dgumbel(trf(f))
          d * (y/p - (1 - y)/(1 - p))
      }, 
      loss = function(y, f){
          trf <- function(f) {
            pmax(pmin(f, 15), -2.5)
          }
          p <- pgumbel(trf(f))
          y <- (y + 1)/2
          -y * log(p) - (1 - y) * log(1 - p)
      }, 
     offset = function(y, w){
         p <- weighted.mean(y > 0, w)
         qgumbel(p)
     }, 
     response = function(f){
         p <- pgumbel(trf(f))
         return(p)
     }, 
  rclass = function(f) (f > 0) + 1, check_y = biny, 
  name = paste("Negative Binomial Likelihood -- Gumbel-Link")))
}


###################################################
### code chunk number 2: IPCW
###################################################

## To account for right-censored observations in the estimation of conditional
## transformation models, we use the inverse probability of censoring weighting
## approach suggested by Graf et al. (1999).

## ngrid: grid of sorted and unique observed survival times
## sT: vector of observed survival times
## event: binary censoring indicator
## G: Kaplan-Meier estimator of the censoring distribution
inv_weights_graf <- function(ngrid, sT, event, G){
  invprob_weights <- c()
  for(i in 1:length(ngrid)){
      r <- ngrid[i] ## grid-point
      weight <- c()
        for(j in 1:length(sT)){ ## cycle through the event/censoring times
          if(event[j] == 1){
             ti <- sT[j] 
             kk <- which(sT[j] == G$time)
             if(length(kk) == 0){ kk <- which.min(abs(sT[j] - G$time))}
             rr <- which(r == G$time)
             ## solve boundary and calculation problems
             if(length(rr) == 0){ rr <- which.min(abs(r - G$time))}
             if(G$surv[kk] == 0) weight_kk <- 60
             else{ weight_kk <- 1/G$surv[kk]}
             if(G$surv[rr] == 0) weight_rr <- 60              
             else{ weight_rr <- 1/G$surv[rr]}
             weight <- c(weight, ifelse(ti <= r, weight_kk, weight_rr))}
                       ## Event time before r => weight = 1/G(ti)
                       ## Event time after r => weight = 1/G(r)
          else{
               ti <- sT[j] ## censoring time
               rr <- which(r == G$time)
               if(length(rr) == 0){ rr <- which.min(abs(r - G$time))}
               if(G$surv[rr] == 0) weight_rr <- 60
               else{ weight_rr <- 1/G$surv[rr]}
               weight <- c(weight, ifelse(ti <= r, 0, weight_rr))}  
                         ## Censoring after r => weight = 1/G(r)
                         ## Censoring before r => weight = 0
           
          }
          weight <- weight*length(sT)/sum(weight) 
          invprob_weights <- c(invprob_weights, weight)
    }
    return(invprob_weights)
}


###################################################
### code chunk number 3: log score
###################################################

## To evaluate the performance of the conditional transformation models and the
## performance of alternative estimation techniques, we use the log score or the
## censored log score, respectively:

## pis: vector of estimated survival probabilites 
log_score <- function(sT, ngrid, pis){
  y <- as.numeric(I(sT <= ngrid)) ## binary survival status 
  sum_it <- 0
  for(i in 1:length(y)){
      ## solve problematic constellations
      if(y[i] == 1 & pis[i] == 0) sum_it <- sum_it + log(1e-8)
      if(y[i] == 0 & pis[i] == 1) sum_it <- sum_it + log(1e-8)
      if(!(y[i] == 1 & pis[i] == 0) & !(y[i] == 0 & pis[i] == 1)){
        if(y[i] == 1) sum_it <- sum_it + log(pis[i])  
        else{ sum_it <- sum_it + log(1 - pis[i])}
      }
  }
  return(sum_it)    
}

log_score_cens <- function(sT, ngrid, pis, weight){
  y <- as.numeric(I(sT <= ngrid)) 
  sum_it <- 0
  for(i in 1:length(y)){      
      if(weight[i] == 0) sum_it <- sum_it
      else{
         ##solve problematic constellations
         if(y[i] == 1 & pis[i] == 0) sum_it <- sum_it + log(1e-8)
         if(y[i] == 0 & pis[i] == 1) sum_it <- sum_it + log(1e-8)
         if(!(y[i] == 1 & pis[i] == 0) & !(y[i] == 0 & pis[i] == 1)){
           if(y[i] == 1) sum_it <- sum_it + log(pis[i])*weight[i]  
           else{ sum_it <- sum_it + log(1 - pis[i])* weight[i]}}
      }
  }
  return(sum_it)    
}


###################################################
### code chunk number 4: sim-prop-hazards
###################################################

## We start with the first simulation example with proportional hazards between
##  the treatment groups. Therefore, we choose Weibull distributed event times 
## for both treatment groups with parameters b_1 = 1 and c_1 = 3 for treatment
## group G1 and b_2 = 1.5 and c_2 = 3 for treatment group G2. The censoring times
## are simulated from a uniform distribution U[0,\gamma], whereby \gamma is
## chosen such that 5\%, 10\%, 25\% and 50\% of censored observations result:

library("ctmDevel")
library("survival")
library("prodlim")

## Simulate 200 observations per treatment group.
n <- 200
## true hazard function Weibull distribution:
h <- function(t, b, c){
     c / (b^c) * t^(c-1)
}
## true survivor function Weibull distribution:
Strue <- function(t, b, c){
   exp(- b^(-c) * t^c)
}


##******************** Simulate data sets ***********************

generate_data <- function(seed, n, gamma1, gamma2, b1, c1, b2, c2){
     set.seed(seed)
     T1 <- rweibull(n, scale = b1, shape = c1)
     C1 <- runif(n, min = 0, max = gamma1)
     T2 <- rweibull(n, scale = b2, shape = c2)
     C2 <- runif(n, min = 0, max = gamma2)
     T <- c(T1, T2) ## event times
     C <- c(C1, C2) ## censoring times
     S <- cbind(T = T, C = C)
     sT <- apply(S, 1, min) ## observed survival times
     df <- data.frame(sT = sT, event = I(T <= C), 
           G = factor(c(rep("G1", n), rep("G2",n)), levels = c("G1", "G2")))
     return(df)
}

## 5% censored observations
df_5 <- generate_data(seed = 12, n = 200, gamma1 = 15, gamma2 = 28, 
                      b1 = 1, c1 = 3, b2 = 1.5, c2 = 3)
## 10% censored observations
df_10 <- generate_data(seed = 12, n = 200, gamma1 = 10, gamma2 = 14, 
                       b1 = 1, c1 = 3, b2 = 1.5, c2 = 3)
## 25% censored observations
df_25 <- generate_data(seed = 12, n = 200, gamma1 = 4, gamma2 = 6, 
                       b1 = 1, c1 = 3, b2 = 1.5, c2 = 3)
## 50% censored observations
df_50 <- generate_data(seed = 12, n = 200, gamma1 = 1.75, gamma2 = 2.6, 
                       b1 = 1, c1 = 3, b2 = 1.5, c2 = 3)


##************ Alternative strategies *********************************

## Calculate Kaplan-Meier estimates
km5 <- prodlim(Surv(sT, event) ~ G, data = df_5)
km10 <- prodlim(Surv(sT, event) ~ G, data = df_10)
km25 <- prodlim(Surv(sT, event) ~ G, data = df_25)
km50 <- prodlim(Surv(sT, event) ~ G, data = df_50)


## Estimate Cox models and calculate corresponding survivor functions
cox_mod5 <- coxph(Surv(sT, event) ~ G, data = df_5)
S_cox5 <- survfit(cox_mod5, newdata = data.frame(sT = df_5$sT, 
                  event = df_5$event, G = df_5$G))
cox_mod10 <- coxph(Surv(sT, event) ~ G, data = df_10)
S_cox10 <- survfit(cox_mod10, newdata = data.frame(sT = df_10$sT, 
                   event = df_10$event, G = df_10$G))
cox_mod25 <- coxph(Surv(sT, event) ~ G, data = df_25)
S_cox25 <- survfit(cox_mod25, newdata = data.frame(sT = df_25$sT, 
                   event = df_25$event, G = df_25$G))
cox_mod50 <- coxph(Surv(sT, event) ~ G, data = df_50)
S_cox50 <- survfit(cox_mod50, newdata = data.frame(sT = df_50$sT, 
                   event = df_50$event, G = df_50$G))



##************** Conditional transformation models ******************

ngrid5 <- sort(unique(df_5$sT))
ngrid10 <- sort(unique(df_10$sT))
ngrid25 <- sort(unique(df_25$sT))
ngrid50 <- sort(unique(df_50$sT))


## Estimate Kaplan-Meier estimates of the censoring distribution
G5 <- prodlim(Surv(sT, event) ~ 1, data = df_5, reverse = TRUE)
G10 <- prodlim(Surv(sT, event) ~ 1, data = df_10, reverse = TRUE)
G25 <- prodlim(Surv(sT, event) ~ 1, data = df_25, reverse = TRUE)
G50 <- prodlim(Surv(sT, event) ~ 1, data = df_50, reverse = TRUE)

## Calculate inverse probability of censoring weights
ipcw5 <- inv_weights_graf(ngrid = ngrid5, sT = df_5$sT, 
                          event = df_5$event, G = G5)
ipcw10 <- inv_weights_graf(ngrid = ngrid10, sT = df_10$sT, 
                           event = df_10$event, G = G10)
ipcw25 <- inv_weights_graf(ngrid = ngrid25, sT = df_25$sT, 
                           event = df_25$event, G = G25)
ipcw50 <- inv_weights_graf(ngrid = ngrid50, sT = df_50$sT, 
                           event = df_50$event, G = G50)

## set up model formula: 
## A separate smooth function over time is estimated for each treatment group;
## constraint = "increasing" guarantees monotone increasing function estimates
fm <- "bbs(sT, df = 4, constraint = \"increasing\") ~ bols(G, df = 2)"
fm <- as.formula(fm)

## Estimate models
ctm_mod5 <- mctm(fm, data = df_5, family = Binomial_gumbel(), 
                 ngrid = ngrid5, weights = ipcw5, control = 
                 boost_control(nu = 0.5, mstop = 3000,
                 trace = TRUE))
ctm_mod10 <- mctm(fm, data = df_10, family = Binomial_gumbel(), 
                  ngrid = ngrid10, weights = ipcw10, control = 
                  boost_control(nu = 0.5, mstop = 3000,
                  trace = TRUE))
ctm_mod25 <- mctm(fm, data = df_25, family = Binomial_gumbel(), 
                  ngrid = ngrid25, weights = ipcw25, control = 
                  boost_control(nu = 0.5, mstop = 3000,
                  trace = TRUE))
ctm_mod50 <- mctm(fm, data = df_50, family = Binomial_gumbel(), 
                  ngrid = ngrid50, weights = ipcw50, control = 
                  boost_control(nu = 0.5, mstop = 3000,
                  trace = TRUE))


## ******************* CTM predictions
## Predict survival probabilites for treatment groups G1 and G2 and 
## arrange them in a matrix with rows corresponding to individuals 
## and columns corresponding to grid points. 
pred_G1_5 <- predict(ctm_mod5, newdata = data.frame(sT = 
                     ctm_mod5$uresponse, G = factor("G1", 
                     levels = c("G1", "G2"))), type = "response")
pred_G1_mat_5 <- matrix(pred_G1_5, ncol = length(ngrid5))
pred_G2_5 <- predict(ctm_mod5, newdata = data.frame(sT = 
                     ctm_mod5$uresponse, G = factor("G2", 
                     levels = c("G1", "G2"))), type = "response")
pred_G2_mat_5 <- matrix(pred_G2_5, ncol = length(ngrid5))

pred_G1_10 <- predict(ctm_mod10, newdata = data.frame(sT = 
                      ctm_mod10$uresponse, G = factor("G1", 
                      levels = c("G1", "G2"))), type = "response")
pred_G1_mat_10 <- matrix(pred_G1_10, ncol = length(ngrid10))
pred_G2_10 <- predict(ctm_mod10, newdata = data.frame(sT = 
                      ctm_mod10$uresponse, G = factor("G2", 
                      levels = c("G1", "G2"))), type = "response")
pred_G2_mat_10 <- matrix(pred_G2_10, ncol = length(ngrid10))

pred_G1_25 <- predict(ctm_mod25, newdata = data.frame(sT = 
                      ctm_mod25$uresponse, G = factor("G1", 
                      levels = c("G1", "G2"))), type = "response")
pred_G1_mat_25 <- matrix(pred_G1_25, ncol = length(ngrid25))
pred_G2_25 <- predict(ctm_mod25, newdata = data.frame(sT = 
                      ctm_mod25$uresponse, G = factor("G2", 
                      levels = c("G1", "G2"))), type = "response")
pred_G2_mat_25 <- matrix(pred_G2_25, ncol = length(ngrid25))

pred_G1_50 <- predict(ctm_mod50, newdata = data.frame(sT = 
                      ctm_mod50$uresponse, G = factor("G1", 
                      levels = c("G1", "G2"))), type = "response")
pred_G1_mat_50 <- matrix(pred_G1_50, ncol = length(ngrid50))
pred_G2_50 <- predict(ctm_mod50, newdata = data.frame(sT = 
                      ctm_mod50$uresponse, G = factor("G2", 
                      levels = c("G1", "G2"))), type = "response")
pred_G2_mat_50 <- matrix(pred_G2_50, ncol = length(ngrid50))


##********************* Evaluation **********************************


## ******************** Calculate MADs

## To get the main paper results the MADs need to be multiplied by 100!
## ctm_G1, ctm_G2: CTM predictions
## S_cox: survivor function from the Cox model
## km: Kaplan-Meier estimator
MADs <- function(ctm_G1, ctm_G2, S_cox, km, ngrid, n, df, 
                 b1, c1, b2, c2){
  MAD_G1_CTM <- sum(abs(ctm_G1[1,] - (1 - Strue(ngrid, b = b1, c = c1))))/
                length(ngrid)
  MAD_G2_CTM <- sum(abs(ctm_G2[1,] - (1 - Strue(ngrid, b = b2, c = c2))))/
                length(ngrid)
  surv_cox <- S_cox$surv
  MAD_G1_cox <- sum(abs(surv_cox[,1] - Strue(ngrid, b = b1, c = c1)))/
                length(ngrid)  
  MAD_G2_cox <- sum(abs(surv_cox[,2*n] - Strue(ngrid, b = b2, c = c2)))/
                length(ngrid)
  pred_km <- predict(km, times = ngrid, newdata = df)
  pred_G1_km <- pred_km[[1]]
  pred_G2_km <- pred_km[[nrow(df)]]
  pred_G1_km[is.na(pred_G1_km)] <- 0
  pred_G2_km[is.na(pred_G2_km)] <- 0
  MAD_G1_km <- sum(abs(pred_G1_km - Strue(ngrid, b = b1, c = c1)))/
               length(ngrid)  
  MAD_G2_km <- sum(abs(pred_G2_km - Strue(ngrid, b = b2, c = c2)))/
               length(ngrid)
  
  MADs <- matrix(c(MAD_G1_CTM, MAD_G2_CTM, MAD_G1_cox, MAD_G2_cox, 
                   MAD_G1_km, MAD_G2_km), nrow = 2, ncol = 3)
  rownames(MADs) <- c("G1", "G2")
  colnames(MADs) <- c("CTM", "Cox", "KM")
  return(MADs)
}

## 5% censored observations
MADs_5 <- MADs(ctm_G1 = pred_G1_mat_5, ctm_G2 = pred_G2_mat_5, 
               S_cox = S_cox5, km = km5, ngrid = ngrid5, 
               n = 200, df = df_5, b1 = 1, c1 = 3, b2 = 1.5, 
               c2 = 3)
## 10% censored observations
MADs_10 <- MADs(ctm_G1 = pred_G1_mat_10, ctm_G2 = pred_G2_mat_10, 
                S_cox = S_cox10, km = km10, ngrid = ngrid10, 
                n = 200, df = df_10, b1 = 1, c1 = 3, b2 = 1.5, 
                c2 = 3)
## 25% censored observations
MADs_25 <- MADs(ctm_G1 = pred_G1_mat_25, ctm_G2 = pred_G2_mat_25, 
                S_cox = S_cox25, km = km25, ngrid = ngrid25, 
                n = 200, df = df_25, b1 = 1, c1 = 3, b2 = 1.5, 
                c2 = 3)
## 50% censored observations
MADs_50 <- MADs(ctm_G1 = pred_G1_mat_50, ctm_G2 = pred_G2_mat_50, 
                S_cox = S_cox50, km = km50, ngrid = ngrid50, 
                n = 200, df = df_50, b1 = 1, c1 = 3, b2 = 1.5, 
                c2 = 3)


## *********************** Calculate uncensored log score
## Simulate new observations 
set.seed(12345) 
sT_G1_test <- rweibull(100, scale = 1, shape = 3)
sT_G2_test <- rweibull(100, scale = 1.5, shape = 3) 

## Log Score CTM G1
log_score_G1_CTM5 <- log_score_G1_CTM10 <- log_score_G1_CTM25 <- log_score_G1_CTM50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_CTM5 <- c(log_score_G1_CTM5, 
                          log_score(sT_G1_test[i], ngrid5,
                                    pred_G1_mat_5[1,]))
   log_score_G1_CTM10 <- c(log_score_G1_CTM10, 
                           log_score(sT_G1_test[i], ngrid10,
                                     pred_G1_mat_10[1,]))
   log_score_G1_CTM25 <- c(log_score_G1_CTM25, 
                           log_score(sT_G1_test[i], ngrid25,
                                     pred_G1_mat_25[1,]))
   log_score_G1_CTM50 <- c(log_score_G1_CTM50, 
                           log_score(sT_G1_test[i], ngrid50,
                                     pred_G1_mat_50[1,]))
}
log_score_G1_CTM5 <- sum(log_score_G1_CTM5)/length(sT_G1_test)
log_score_G1_CTM10 <- sum(log_score_G1_CTM10)/length(sT_G1_test)
log_score_G1_CTM25 <- sum(log_score_G1_CTM25)/length(sT_G1_test)
log_score_G1_CTM50 <- sum(log_score_G1_CTM50)/length(sT_G1_test)

## Log Score CTM G2
log_score_G2_CTM5 <- log_score_G2_CTM10 <- log_score_G2_CTM25 <- log_score_G2_CTM50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_CTM5 <- c(log_score_G2_CTM5, 
                          log_score(sT_G2_test[i], ngrid5,
                                    pred_G2_mat_5[1,]))
   log_score_G2_CTM10 <- c(log_score_G2_CTM10, 
                           log_score(sT_G2_test[i], ngrid10,
                                     pred_G2_mat_10[1,]))
   log_score_G2_CTM25 <- c(log_score_G2_CTM25, 
                           log_score(sT_G2_test[i], ngrid25,
                                     pred_G2_mat_25[1,]))
   log_score_G2_CTM50 <- c(log_score_G2_CTM50, 
                           log_score(sT_G2_test[i], ngrid50,
                                     pred_G2_mat_50[1,]))
}
log_score_G2_CTM5 <- sum(log_score_G2_CTM5)/length(sT_G2_test)
log_score_G2_CTM10 <- sum(log_score_G2_CTM10)/length(sT_G2_test)
log_score_G2_CTM25 <- sum(log_score_G2_CTM25)/length(sT_G2_test)
log_score_G2_CTM50 <- sum(log_score_G2_CTM50)/length(sT_G2_test)

## Log Score Cox G1
log_score_G1_cox5 <- log_score_G1_cox10 <- log_score_G1_cox25 <- log_score_G1_cox50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_cox5 <- c(log_score_G1_cox5, 
                          log_score(sT_G1_test[i], ngrid5,
                                    1 - S_cox5$surv[,1]))
   log_score_G1_cox10 <- c(log_score_G1_cox10, 
                           log_score(sT_G1_test[i], ngrid10,
                                     1 - S_cox10$surv[,1]))
   log_score_G1_cox25 <- c(log_score_G1_cox25, 
                           log_score(sT_G1_test[i], ngrid25,
                                     1 - S_cox25$surv[,1]))
   log_score_G1_cox50 <- c(log_score_G1_cox50, 
                           log_score(sT_G1_test[i], ngrid50,
                                     1 - S_cox50$surv[,1]))
}
log_score_G1_cox5 <- sum(log_score_G1_cox5)/length(sT_G1_test)
log_score_G1_cox10 <- sum(log_score_G1_cox10)/length(sT_G1_test)
log_score_G1_cox25 <- sum(log_score_G1_cox25)/length(sT_G1_test)
log_score_G1_cox50 <- sum(log_score_G1_cox50)/length(sT_G1_test)

## Log Score Cox G2  
log_score_G2_cox5 <- log_score_G2_cox10 <- log_score_G2_cox25 <- log_score_G2_cox50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_cox5 <- c(log_score_G2_cox5, 
                          log_score(sT_G2_test[i], ngrid5,
                                    1 - S_cox5$surv[,2*n]))
   log_score_G2_cox10 <- c(log_score_G2_cox10, 
                           log_score(sT_G2_test[i], ngrid10,
                                     1 - S_cox10$surv[,2*n]))
   log_score_G2_cox25 <- c(log_score_G2_cox25, 
                           log_score(sT_G2_test[i], ngrid25,
                                     1 - S_cox25$surv[,2*n]))
   log_score_G2_cox50 <- c(log_score_G2_cox50, 
                           log_score(sT_G2_test[i], ngrid50,
                                     1 - S_cox50$surv[,2*n]))
}
log_score_G2_cox5 <- sum(log_score_G2_cox5)/length(sT_G2_test)
log_score_G2_cox10 <- sum(log_score_G2_cox10)/length(sT_G2_test)
log_score_G2_cox25 <- sum(log_score_G2_cox25)/length(sT_G2_test)
log_score_G2_cox50 <- sum(log_score_G2_cox50)/length(sT_G2_test)

## Log Score KM G1
pred_km5 <- predict(km5, times = ngrid5, newdata = df_5)
pred_G1_km5 <- pred_km5[[1]]
pred_G1_km5[is.na(pred_G1_km5)] <- 0
pred_km10 <- predict(km10, times = ngrid10, newdata = df_10)
pred_G1_km10 <- pred_km10[[1]]
pred_G1_km10[is.na(pred_G1_km10)] <- 0
pred_km25 <- predict(km25, times = ngrid25, newdata = df_25)
pred_G1_km25 <- pred_km25[[1]]
pred_G1_km25[is.na(pred_G1_km25)] <- 0
pred_km50 <- predict(km50, times = ngrid50, newdata = df_50)
pred_G1_km50 <- pred_km50[[1]]
pred_G1_km50[is.na(pred_G1_km50)] <- 0

log_score_G1_km5 <- log_score_G1_km10 <- log_score_G1_km25 <- 
log_score_G1_km50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_km5 <- c(log_score_G1_km5, 
                         log_score(sT_G1_test[i], ngrid5,
                                   1 - pred_G1_km5))
   log_score_G1_km10 <- c(log_score_G1_km10, 
                          log_score(sT_G1_test[i], ngrid10,
                                    1 - pred_G1_km10))
   log_score_G1_km25 <- c(log_score_G1_km25, 
                          log_score(sT_G1_test[i], ngrid25,
                                    1 - pred_G1_km25))
   log_score_G1_km50 <- c(log_score_G1_km50, 
                          log_score(sT_G1_test[i], ngrid50,
                                    1 - pred_G1_km50))
}
log_score_G1_km5 <- sum(log_score_G1_km5)/length(sT_G1_test)
log_score_G1_km10 <- sum(log_score_G1_km10)/length(sT_G1_test)
log_score_G1_km25 <- sum(log_score_G1_km25)/length(sT_G1_test)
log_score_G1_km50 <- sum(log_score_G1_km50)/length(sT_G1_test)


## Log Score KM G2 
pred_G2_km5 <- pred_km5[[nrow(df_5)]]
pred_G2_km5[is.na(pred_G2_km5)] <- 0
pred_G2_km10 <- pred_km10[[nrow(df_10)]]
pred_G2_km10[is.na(pred_G2_km10)] <- 0
pred_G2_km25 <- pred_km25[[nrow(df_25)]]
pred_G2_km25[is.na(pred_G2_km25)] <- 0
pred_G2_km50 <- pred_km50[[nrow(df_50)]]
pred_G2_km50[is.na(pred_G2_km50)] <- 0 
log_score_G2_km5 <- log_score_G2_km10 <- log_score_G2_km25 <- 
log_score_G2_km50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_km5 <- c(log_score_G2_km5, 
                         log_score(sT_G2_test[i], ngrid5,
                                   1 - pred_G2_km5))
   log_score_G2_km10 <- c(log_score_G2_km10, 
                          log_score(sT_G2_test[i], ngrid10,
                                    1 - pred_G2_km10))
   log_score_G2_km25 <- c(log_score_G2_km25, 
                          log_score(sT_G2_test[i], ngrid25,
                                    1 - pred_G2_km25))
   log_score_G2_km50 <- c(log_score_G2_km50, 
                          log_score(sT_G2_test[i], ngrid50,
                                    1 - pred_G2_km50))
}
log_score_G2_km5 <- sum(log_score_G2_km5)/length(sT_G2_test)
log_score_G2_km10 <- sum(log_score_G2_km10)/length(sT_G2_test)
log_score_G2_km25 <- sum(log_score_G2_km25)/length(sT_G2_test)
log_score_G2_km50 <- sum(log_score_G2_km50)/length(sT_G2_test)


###################################################
### code chunk number 5: sim-nonprop-hazards
###################################################

## In the second simulation example we assume non-proportional hazards between
## the treatment groups. Therefore, we choose Weibull distributed event times for
## both treatment groups with parameters b_1 = 1.5 and c_1 = 3 for treatment
## group G1 and b_2 = 1 and c_2 = 1 for treatment group G2. The censoring times
## are simulated from a uniform distribution U[0,\gamma], whereby \gamma is
## chosen such that 5\%, 10\%, 25\% and 50\% of censored observations result:


##*************** Simulate data sets ******************
## 5% censored observations
df_5 <- generate_data(seed = 12, n = 200, gamma1 = 22, gamma2 = 16, 
                      b1 = 1.5, c1 = 3, b2 = 1, c2 = 1)
## 10% censored observations
df_10 <- generate_data(seed = 12, n = 200, gamma1 = 14, gamma2 = 8.5, 
                       b1 = 1.5, c1 = 3, b2 = 1, c2 = 1)
## 25% censored observations
df_25 <- generate_data(seed = 12, n = 200, gamma1 = 6, gamma2 = 4.3, 
                       b1 = 1.5, c1 = 3, b2 = 1, c2 = 1)
## 50% censored observations
df_50 <- generate_data(seed = 12, n = 200, gamma1 = 2.6, gamma2 = 1.65, 
                       b1 = 1.5, c1 = 3, b2 = 1, c2 = 1)


##******************* Alternative Strategies ***************

## Kaplan-Meier estimators
km5 <- prodlim(Surv(sT, event) ~ G, data = df_5)
km10 <- prodlim(Surv(sT, event) ~ G, data = df_10)
km25 <- prodlim(Surv(sT, event) ~ G, data = df_25)
km50 <- prodlim(Surv(sT, event) ~ G, data = df_50)

## Cox Models and corresponding survivor functions
cox_mod5 <- coxph(Surv(sT, event) ~ G, data = df_5)
S_cox5 <- survfit(cox_mod5, newdata = data.frame(sT = df_5$sT, 
                            event = df_5$event, G = df_5$G))
cox_mod10 <- coxph(Surv(sT, event) ~ G, data = df_10)
S_cox10 <- survfit(cox_mod10, newdata = data.frame(sT = df_10$sT, 
                              event = df_10$event, G = df_10$G))
cox_mod25 <- coxph(Surv(sT, event) ~ G, data = df_25)
S_cox25 <- survfit(cox_mod25, newdata = data.frame(sT = df_25$sT, 
                              event = df_25$event, G = df_25$G))
cox_mod50 <- coxph(Surv(sT, event) ~ G, data = df_50)
S_cox50 <- survfit(cox_mod50, newdata = data.frame(sT = df_50$sT, 
                              event = df_50$event, G = df_50$G))
## Stratified Cox-Model and corresponding survivor functions
cox_mod5.2 <- coxph(Surv(sT, event) ~ G + strata(G), data = df_5)
S_cox5.2 <- survfit(cox_mod5.2, newdata = data.frame(sT = df_5$sT, 
                                event = df_5$event, G = df_5$G))
cox_mod10.2 <- coxph(Surv(sT, event) ~ G + strata(G), data = df_10)
S_cox10.2 <- survfit(cox_mod10.2, newdata = data.frame(sT = df_10$sT, 
                                  event = df_10$event, G = df_10$G))
cox_mod25.2 <- coxph(Surv(sT, event) ~ G + strata(G), data = df_25)
S_cox25.2 <- survfit(cox_mod25.2, newdata = data.frame(sT = df_25$sT, 
                                  event = df_25$event, G = df_25$G))
cox_mod50.2 <- coxph(Surv(sT, event) ~ G + strata(G), data = df_50)
S_cox50.2 <- survfit(cox_mod50.2, newdata = data.frame(sT = df_50$sT, 
                                  event = df_50$event, G = df_50$G))


##************** Conditional transformation models

ngrid5 <- sort(unique(df_5$sT))
ngrid10 <- sort(unique(df_10$sT))
ngrid25 <- sort(unique(df_25$sT))
ngrid50 <- sort(unique(df_50$sT))

## Kaplan-Meier estimators of the censoring distribution
G5 <- prodlim(Surv(sT, event) ~ 1, data = df_5, reverse = TRUE)
G10 <- prodlim(Surv(sT, event) ~ 1, data = df_10, reverse = TRUE)
G25 <- prodlim(Surv(sT, event) ~ 1, data = df_25, reverse = TRUE)
G50 <- prodlim(Surv(sT, event) ~ 1, data = df_50, reverse = TRUE)

## Inverse probability of censoring weights
ipcw5 <- inv_weights_graf(ngrid = ngrid5, sT = df_5$sT, 
                          event = df_5$event, G = G5)
ipcw10 <- inv_weights_graf(ngrid = ngrid10, sT = df_10$sT,
                           event = df_10$event, G = G10)
ipcw25 <- inv_weights_graf(ngrid = ngrid25, sT = df_25$sT,
                           event = df_25$event, G = G25)
ipcw50 <- inv_weights_graf(ngrid = ngrid50, sT = df_50$sT,
                           event = df_50$event, G = G50)
## Set up model formula
fm <- "bbs(sT, df = 8, constraint = \"increasing\") ~ bols(G, df = 2)"
fm <- as.formula(fm)

## Estimate models
ctm_mod5 <- mctm(fm, data = df_5, family = Binomial_gumbel(), 
                 ngrid = ngrid5, weights = ipcw5, control = 
                 boost_control(nu = 0.8, mstop = 3000,
                 trace = TRUE))
ctm_mod10 <- mctm(fm, data = df_10, family = Binomial_gumbel(), 
                  ngrid = ngrid10, weights = ipcw10, control = 
                  boost_control(nu = 0.5, mstop = 3000,
                  trace = TRUE))
ctm_mod25 <- mctm(fm, data = df_25, family = Binomial_gumbel(), 
                  ngrid = ngrid25, weights = ipcw25, control = 
                  boost_control(nu = 0.5, mstop = 3000,
                  trace = TRUE))
ctm_mod50 <- mctm(fm, data = df_50, family = Binomial_gumbel(), 
                  ngrid = ngrid50, weights = ipcw50, control = 
                  boost_control(nu = 0.5, mstop = 3000,
                  trace = TRUE))


##*************** CTM predictions

pred_G1_5 <- predict(ctm_mod5, newdata = data.frame(sT = 
                     ctm_mod5$uresponse, G = factor("G1", 
                     levels = c("G1", "G2"))), type = "response")
pred_G1_mat_5 <- matrix(pred_G1_5, ncol = length(ngrid5))
pred_G2_5 <- predict(ctm_mod5, newdata = data.frame(sT = 
                     ctm_mod5$uresponse, G = factor("G2", 
                     levels = c("G1", "G2"))), type = "response")
pred_G2_mat_5 <- matrix(pred_G2_5, ncol = length(ngrid5))

pred_G1_10 <- predict(ctm_mod10, newdata = data.frame(sT = 
                      ctm_mod10$uresponse, G = factor("G1", 
                      levels = c("G1", "G2"))), type = "response")
pred_G1_mat_10 <- matrix(pred_G1_10, ncol = length(ngrid10))
pred_G2_10 <- predict(ctm_mod10, newdata = data.frame(sT = 
                      ctm_mod10$uresponse, G = factor("G2", 
                      levels = c("G1", "G2"))), type = "response")
pred_G2_mat_10 <- matrix(pred_G2_10, ncol = length(ngrid10))

pred_G1_25 <- predict(ctm_mod25, newdata = data.frame(sT = 
                      ctm_mod25$uresponse, G = factor("G1", 
                      levels = c("G1", "G2"))), type = "response")
pred_G1_mat_25 <- matrix(pred_G1_25, ncol = length(ngrid25))
pred_G2_25 <- predict(ctm_mod25, newdata = data.frame(sT = 
                      ctm_mod25$uresponse, G = factor("G2", 
                      levels = c("G1", "G2"))), type = "response")
pred_G2_mat_25 <- matrix(pred_G2_25, ncol = length(ngrid25))

pred_G1_50 <- predict(ctm_mod50, newdata = data.frame(sT = 
                      ctm_mod50$uresponse, G = factor("G1", 
                      levels = c("G1", "G2"))), type = "response")
pred_G1_mat_50 <- matrix(pred_G1_50, ncol = length(ngrid50))
pred_G2_50 <- predict(ctm_mod50, newdata = data.frame(sT = 
                      ctm_mod50$uresponse, G = factor("G2", 
                      levels = c("G1", "G2"))), type = "response")
pred_G2_mat_50 <- matrix(pred_G2_50, ncol = length(ngrid50))
 


##************************ Evaluation *********************************

##********* Calculate MADs
## To get the main paper results the MADs need to be multiplied by 100!

## 5% censored observations
MADs_5 <- MADs(ctm_G1 = pred_G1_mat_5, ctm_G2 = pred_G2_mat_5, 
               S_cox = S_cox5, km = km5, ngrid = ngrid5, 
               n = 200, df = df_5, b1 = 1.5, c1 = 3, b2 = 1, 
               c2 = 1)
## 10% censored observations
MADs_10 <- MADs(ctm_G1 = pred_G1_mat_10, ctm_G2 = pred_G2_mat_10, 
                S_cox = S_cox10, km = km10, ngrid = ngrid10, 
                n = 200, df = df_10, b1 = 1.5, c1 = 3, b2 = 1, 
                c2 = 1)
## 25% censored observations
MADs_25 <- MADs(ctm_G1 = pred_G1_mat_25, ctm_G2 = pred_G2_mat_25, 
                S_cox = S_cox25, km = km25, ngrid = ngrid25, 
                n = 200, df = df_25, b1 = 1.5, c1 = 3, b2 = 1, 
                c2 = 1)
## 50% censored observations
MADs_50 <- MADs(ctm_G1 = pred_G1_mat_50, ctm_G2 = pred_G2_mat_50, 
                S_cox = S_cox50, km = km50, ngrid = ngrid50, 
                n = 200, df = df_50, b1 = 1.5, c1 = 3, b2 = 1, 
                c2 = 1)

## MAD Stratified Cox
## Calculation is more complicated since the survfit-function does not 
## return the desired result.
MAD_strat_cox <- function(strat_cox_mod, ngrid, b1, c1, b2, c2){
  baseline <- basehaz(strat_cox_mod, centered = FALSE)
  baselineG1 <- baseline[baseline$strata == "G=G1",]
  baselineG2 <- baseline[baseline$strata == "G=G2",]
  survG1_strata <- survG2_strata <- rep(0, length(ngrid))
  for(i in 1:length(ngrid)){
    k <- max(which(baselineG1$time <= ngrid[i]))
    l <- max(which(baselineG2$time <= ngrid[i]))
    ## Do not worry about the warnings. The next lines solve 
    ## the problem!
    if(k == -Inf) k <- 1
    if(l == -Inf) l <- 1
    survG1_strata[i] <- exp(-baselineG1$hazard[k])
    survG2_strata[i] <- exp(-baselineG2$hazard[l])
  }
  MAD_StratCox_G1 <- sum(abs(survG1_strata - Strue(ngrid, b = b1, 
                         c = c1)))/length(ngrid)   
  MAD_StratCox_G2 <- sum(abs(survG2_strata - Strue(ngrid, b = b2, 
                         c = c2)))/length(ngrid)
  return(c(MAD_StratCox_G1, MAD_StratCox_G2))
}
MAD_StratCox_5 <- MAD_strat_cox(cox_mod5.2, ngrid5, b1 = 1.5, 
                                c1 = 3, b2 = 1, c2 = 1)
MAD_StratCox_10 <- MAD_strat_cox(cox_mod10.2, ngrid10, b1 = 1.5, 
                                 c1 = 3, b2 = 1, c2 = 1)
MAD_StratCox_25 <- MAD_strat_cox(cox_mod25.2, ngrid25, b1 = 1.5, 
                                 c1 = 3, b2 = 1, c2 = 1)
MAD_StratCox_50 <- MAD_strat_cox(cox_mod50.2, ngrid50, b1 = 1.5, 
                                 c1 = 3, b2 = 1, c2 = 1)
  

## ************* Calculate uncensored log score for new observations
set.seed(12345) 
sT_G1_test <- rweibull(100, scale = 1.5, shape = 3)
sT_G2_test <- rweibull(100, scale = 1, shape = 1) 

## Log score CTM G1:
log_score_G1_CTM5 <- log_score_G1_CTM10 <- log_score_G1_CTM25 <-
log_score_G1_CTM50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_CTM5 <- c(log_score_G1_CTM5, 
                          log_score(sT_G1_test[i], ngrid5, 
                                    pred_G1_mat_5[1,]))
   log_score_G1_CTM10 <- c(log_score_G1_CTM10, 
                           log_score(sT_G1_test[i], ngrid10, 
                                     pred_G1_mat_10[1,]))
   log_score_G1_CTM25 <- c(log_score_G1_CTM25, 
                           log_score(sT_G1_test[i], ngrid25, 
                                     pred_G1_mat_25[1,]))
   log_score_G1_CTM50 <- c(log_score_G1_CTM50, 
                           log_score(sT_G1_test[i], ngrid50, 
                                     pred_G1_mat_50[1,]))
}
log_score_G1_CTM5 <- sum(log_score_G1_CTM5)/length(sT_G1_test)
log_score_G1_CTM10 <- sum(log_score_G1_CTM10)/length(sT_G1_test)
log_score_G1_CTM25 <- sum(log_score_G1_CTM25)/length(sT_G1_test)
log_score_G1_CTM50 <- sum(log_score_G1_CTM50)/length(sT_G1_test)

## Log score CTM G2
log_score_G2_CTM5 <- log_score_G2_CTM10 <- log_score_G2_CTM25 <-
log_score_G2_CTM50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_CTM5 <- c(log_score_G2_CTM5, 
                          log_score(sT_G2_test[i], ngrid5, 
                                    pred_G2_mat_5[1,]))
   log_score_G2_CTM10 <- c(log_score_G2_CTM10, 
                           log_score(sT_G2_test[i], ngrid10, 
                                     pred_G2_mat_10[1,]))
   log_score_G2_CTM25 <- c(log_score_G2_CTM25, 
                           log_score(sT_G2_test[i], ngrid25, 
                                     pred_G2_mat_25[1,]))
   log_score_G2_CTM50 <- c(log_score_G2_CTM50, 
                           log_score(sT_G2_test[i], ngrid50, 
                                     pred_G2_mat_50[1,]))
}
log_score_G2_CTM5 <- sum(log_score_G2_CTM5)/length(sT_G2_test)
log_score_G2_CTM10 <- sum(log_score_G2_CTM10)/length(sT_G2_test)
log_score_G2_CTM25 <- sum(log_score_G2_CTM25)/length(sT_G2_test)
log_score_G2_CTM50 <- sum(log_score_G2_CTM50)/length(sT_G2_test)

## Log score Cox G1
log_score_G1_cox5 <- log_score_G1_cox10 <- log_score_G1_cox25 <-
log_score_G1_cox50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_cox5 <- c(log_score_G1_cox5, 
                          log_score(sT_G1_test[i], ngrid5, 
                                    1 - S_cox5$surv[,1]))
   log_score_G1_cox10 <- c(log_score_G1_cox10, 
                           log_score(sT_G1_test[i], ngrid10,
                                     1 - S_cox10$surv[,1]))
   log_score_G1_cox25 <- c(log_score_G1_cox25,
                           log_score(sT_G1_test[i], ngrid25,
                                     1 - S_cox25$surv[,1]))
   log_score_G1_cox50 <- c(log_score_G1_cox50, 
                           log_score(sT_G1_test[i], ngrid50,
                                     1 - S_cox50$surv[,1]))
}
log_score_G1_cox5 <- sum(log_score_G1_cox5)/length(sT_G1_test)
log_score_G1_cox10 <- sum(log_score_G1_cox10)/length(sT_G1_test)
log_score_G1_cox25 <- sum(log_score_G1_cox25)/length(sT_G1_test)
log_score_G1_cox50 <- sum(log_score_G1_cox50)/length(sT_G1_test)

## Log score Cox G2  
log_score_G2_cox5 <- log_score_G2_cox10 <- log_score_G2_cox25 <-
log_score_G2_cox50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_cox5 <- c(log_score_G2_cox5, 
                          log_score(sT_G2_test[i], ngrid5, 
                                    1 - S_cox5$surv[,2*n]))
   log_score_G2_cox10 <- c(log_score_G2_cox10, 
                           log_score(sT_G2_test[i], ngrid10, 
                                     1 - S_cox10$surv[,2*n]))
   log_score_G2_cox25 <- c(log_score_G2_cox25, 
                           log_score(sT_G2_test[i], ngrid25, 
                                     1 - S_cox25$surv[,2*n]))
   log_score_G2_cox50 <- c(log_score_G2_cox50, 
                           log_score(sT_G2_test[i], ngrid50, 
                                     1 - S_cox50$surv[,2*n]))
}
log_score_G2_cox5 <- sum(log_score_G2_cox5)/length(sT_G2_test)
log_score_G2_cox10 <- sum(log_score_G2_cox10)/length(sT_G2_test)
log_score_G2_cox25 <- sum(log_score_G2_cox25)/length(sT_G2_test)
log_score_G2_cox50 <- sum(log_score_G2_cox50)/length(sT_G2_test)

## Log score stratified Cox 
## Calculate survival probabilites for stratified Cox models
surv_probs <- function(strat_cox_mod, ngrid){
  baseline <- basehaz(strat_cox_mod, centered = FALSE)
  baselineG1 <- baseline[baseline$strata == "G=G1",]
  baselineG2 <- baseline[baseline$strata == "G=G2",]
  survG1_strata <- survG2_strata <- rep(0, length(ngrid))
  for(i in 1:length(ngrid)){
    k <- max(which(baselineG1$time <= ngrid[i]))
    l <- max(which(baselineG2$time <= ngrid[i]))
    ## Do not worry about the warnings. The next two lines 
    ## solve the problem!
    if(k == -Inf) k <- 1
    if(l == -Inf) l <- 1
    survG1_strata[i] <- exp(-baselineG1$hazard[k])
    survG2_strata[i] <- exp(-baselineG2$hazard[l])
  }
  surv_probs <- vector(mode = "list", length = 2)
  names(surv_probs) <- c("G1", "G2")
  surv_probs[[1]] <- survG1_strata
  surv_probs[[2]] <- survG2_strata
  return(surv_probs)
}
surv_probs_5 <- surv_probs(cox_mod5.2, ngrid5)
surv_probs_10 <- surv_probs(cox_mod10.2, ngrid10)
surv_probs_25 <- surv_probs(cox_mod25.2, ngrid25)
surv_probs_50 <- surv_probs(cox_mod50.2, ngrid50)

## Log score stratified Cox G1
log_score_G1_stratcox5 <- log_score_G1_stratcox10 <- c()
log_score_G1_stratcox25 <- log_score_G1_stratcox50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_stratcox5 <- c(log_score_G1_stratcox5, 
                               log_score(sT_G1_test[i], ngrid5, 
                                         1 - surv_probs_5[[1]]))
   log_score_G1_stratcox10 <- c(log_score_G1_stratcox10, 
                                log_score(sT_G1_test[i], ngrid10, 
                                          1 - surv_probs_10[[1]]))
   log_score_G1_stratcox25 <- c(log_score_G1_stratcox25, 
                                log_score(sT_G1_test[i], ngrid25, 
                                          1 - surv_probs_25[[1]]))
   log_score_G1_stratcox50 <- c(log_score_G1_stratcox50, 
                                log_score(sT_G1_test[i], ngrid50, 
                                          1 - surv_probs_50[[1]]))
}
log_score_G1_stratcox5 <- sum(log_score_G1_stratcox5)/length(sT_G1_test)
log_score_G1_stratcox10 <- sum(log_score_G1_stratcox10)/length(sT_G1_test)
log_score_G1_stratcox25 <- sum(log_score_G1_stratcox25)/length(sT_G1_test)
log_score_G1_stratcox50 <- sum(log_score_G1_stratcox50)/length(sT_G1_test)

## Log-Score stratified Cox G2  
log_score_G2_stratcox5 <- log_score_G2_stratcox10 <- c()
log_score_G2_stratcox25 <- log_score_G2_stratcox50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_stratcox5 <- c(log_score_G2_stratcox5, 
                               log_score(sT_G2_test[i], ngrid5, 
                                         1 - surv_probs_5[[2]]))
   log_score_G2_stratcox10 <- c(log_score_G2_stratcox10, 
                                log_score(sT_G2_test[i], ngrid10, 
                                          1 - surv_probs_10[[2]]))
   log_score_G2_stratcox25 <- c(log_score_G2_stratcox25, 
                                log_score(sT_G2_test[i], ngrid25, 
                                          1 - surv_probs_25[[2]]))
   log_score_G2_stratcox50 <- c(log_score_G2_stratcox50, 
                                log_score(sT_G2_test[i], ngrid50, 
                                          1 - surv_probs_50[[2]]))
}
log_score_G2_stratcox5 <- sum(log_score_G2_stratcox5)/length(sT_G2_test)
log_score_G2_stratcox10 <- sum(log_score_G2_stratcox10)/length(sT_G2_test)
log_score_G2_stratcox25 <- sum(log_score_G2_stratcox25)/length(sT_G2_test)
log_score_G2_stratcox50 <- sum(log_score_G2_stratcox50)/length(sT_G2_test)

## Log score KM G1
pred_km5 <- predict(km5, times = ngrid5, newdata = df_5)
pred_G1_km5 <- pred_km5[[1]]
pred_G1_km5[is.na(pred_G1_km5)] <- 0
pred_km10 <- predict(km10, times = ngrid10, newdata = df_10)
pred_G1_km10 <- pred_km10[[1]]
pred_G1_km10[is.na(pred_G1_km10)] <- 0
pred_km25 <- predict(km25, times = ngrid25, newdata = df_25)
pred_G1_km25 <- pred_km25[[1]]
pred_G1_km25[is.na(pred_G1_km25)] <- 0
pred_km50 <- predict(km50, times = ngrid50, newdata = df_50)
pred_G1_km50 <- pred_km50[[1]]
pred_G1_km50[is.na(pred_G1_km50)] <- 0

log_score_G1_km5 <- log_score_G1_km10 <- log_score_G1_km25 <-
log_score_G1_km50 <- c()
for(i in 1:length(sT_G1_test)){
   log_score_G1_km5 <- c(log_score_G1_km5, 
                         log_score(sT_G1_test[i], ngrid5,
                                   1 - pred_G1_km5))
   log_score_G1_km10 <- c(log_score_G1_km10, 
                          log_score(sT_G1_test[i], ngrid10,
                                    1 - pred_G1_km10))
   log_score_G1_km25 <- c(log_score_G1_km25, 
                          log_score(sT_G1_test[i], ngrid25,
                                    1 - pred_G1_km25))
   log_score_G1_km50 <- c(log_score_G1_km50, 
                          log_score(sT_G1_test[i], ngrid50,
                                    1 - pred_G1_km50))
}
log_score_G1_km5 <- sum(log_score_G1_km5)/length(sT_G1_test)
log_score_G1_km10 <- sum(log_score_G1_km10)/length(sT_G1_test)
log_score_G1_km25 <- sum(log_score_G1_km25)/length(sT_G1_test)
log_score_G1_km50 <- sum(log_score_G1_km50)/length(sT_G1_test)

## Log score KM G2
pred_G2_km5 <- pred_km5[[nrow(df_5)]]
pred_G2_km5[is.na(pred_G2_km5)] <- 0
pred_G2_km10 <- pred_km10[[nrow(df_10)]]
pred_G2_km10[is.na(pred_G2_km10)] <- 0
pred_G2_km25 <- pred_km25[[nrow(df_25)]]
pred_G2_km25[is.na(pred_G2_km25)] <- 0
pred_G2_km50 <- pred_km50[[nrow(df_50)]]
pred_G2_km50[is.na(pred_G2_km50)] <- 0 

log_score_G2_km5 <- log_score_G2_km10 <- log_score_G2_km25 <-
log_score_G2_km50 <- c()
for(i in 1:length(sT_G2_test)){
   log_score_G2_km5 <- c(log_score_G2_km5, 
                         log_score(sT_G2_test[i], ngrid5,
                                   1 - pred_G2_km5))
   log_score_G2_km10 <- c(log_score_G2_km10, 
                          log_score(sT_G2_test[i], ngrid10,
                                    1 - pred_G2_km10))
   log_score_G2_km25 <- c(log_score_G2_km25, 
                          log_score(sT_G2_test[i], ngrid25,
                                    1 - pred_G2_km25))
   log_score_G2_km50 <- c(log_score_G2_km50, 
                          log_score(sT_G2_test[i], ngrid50,
                                    1 - pred_G2_km50))
}
log_score_G2_km5 <- sum(log_score_G2_km5)/length(sT_G2_test)
log_score_G2_km10 <- sum(log_score_G2_km10)/length(sT_G2_test)
log_score_G2_km25 <- sum(log_score_G2_km25)/length(sT_G2_test)
log_score_G2_km50 <- sum(log_score_G2_km50)/length(sT_G2_test)

