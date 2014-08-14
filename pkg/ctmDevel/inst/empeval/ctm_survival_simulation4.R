
##********************************************************************
## CTM for survival data: Simulation 4 
##********************************************************************

## The R code given below enables the reader to reproduce the simulation results 
## of Simulation 4 given in Section 3 of the main paper.
## All analysis were carried out in R version 3.1.1. Furthermore, one needs 
## to install the add-on packages ctmDevel (depending on mboostDevel), 
## survival (for the estimation of Cox models), prodlim (for the calculation of
## Kaplan-Meier estimators), and party (for the estimation of conditional
## random forests).

## Load family
source("Binomial_extreme.R")
## The function Binomial_extreme defines a Family-object which is  compatible with
## the families from the mboost-package. The Binomial_extreme-family is used for
## binary response variables and the corresponding link function is set to the
## distribution function of the minimum-extreme value distribution.



## **************************** IPCW
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


## Uncensored integrated log score
## sT: observed survival time.
## ngrid: grid points over which survival probabilities are estimated.
## pis: estimated survival probabilites over the grid points. 
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
  return((-1)*sum_it/length(ngrid))  ## *(-1) since negative Binomial log-likelihood
                                     ## is used in the log score;
                                     ## /length(ngrid) to get log score per observation
                                     ## AND grid point  
}

## load libraries
library("ctmDevel")
library("survival")
library("prodlim")
library("party")

## true hazard function Weibull distribution:
h <- function(t, b, c){
     c / (b^c) * t^(c-1)
}
## true survivor function Weibull distribution:
Strue <- function(t, b, c){
   exp(- b^(-c) * t^c)
}

## ************************* Simulation 4 ********************************************
## Non-proportional hazards
## Influential treatment group G
## Influential continuous covariate x


##******************** Simulate data sets ***********************

## generate_data(): Function to generate data sets.
## seed: set seed.
## n: observations per treatment group.
## p: proportion of censoring, p \in {0.05, 0.1, 0.25, 0.5}.
## b: fixed scale parameter for the Weibull distribution.
## x: grid of x-values.
generate_data <- function(seed, n, p, b, x){
     set.seed(seed)
 
     T1 <- T2 <- c()
     for(i in 1:(length(x)/2)){
       ## c1 = 2 + x_i^2; c2 = 2.5 + x_i^2
       ## = Weibull shape parameters for the treatment groups.
       T1 <- c(T1, rweibull(1, scale = b, shape = (2 + x[i]^2)))
       T2 <- c(T2, rweibull(1, scale = b, shape = (2.5 + x[i]^2)))
     }
     n_cens <- round(n*p, digits = 0)
 
     idx_cens1 <- sample(1:n, n_cens, replace = FALSE)
     C_sub <- sapply(T1[idx_cens1], FUN = function(y){runif(1, min = 0, max = y)})
     C1 <- rep(5000,n)
     C1[idx_cens1] <- C_sub
 
     idx_cens2 <- sample(1:n, n_cens, replace = FALSE)
     C_sub <- sapply(T2[idx_cens2], FUN = function(y){runif(1, min = 0, max = y)})
     C2 <- rep(5000,n)
     C2[idx_cens2] <- C_sub
 
     T <- c(T1, T2)
     C <- c(C1, C2)
     S <- cbind(T = T, C = C)
     sT <- apply(S, 1, min)
     
     df <- data.frame(sT = sT, event = I(T <= C), 
                      G = factor(c(rep("G1", n), rep("G2",n)), levels = c("G1", "G2")),
                      xvar = x)
     return(df)
}


## Choose fixed grid of x-values:
## x ~ U[0,2]
n <- 300 ## Observations per treatment group
x_var <- seq(from = 0, to = 2, length = n)
x_ges <- expand.grid(x = x_var, G = factor(c("G1", "G2"), levels = c("G1", "G2")))

## seeds
set.seed(80)
seeds <- round(runif(100), digits = 4)*10000

## Generate list of data sets
## Fixed Weibull scale parameter b = exp(0.5).
list_df_5 <- list_df_10 <- list_df_25 <- list_df_50 <- vector(mode = "list", 
                                                              length = 100)
for(i in 1:100){
  list_df_5[[i]] <- generate_data(seeds[i], n = n, p = 0.05, b = exp(0.5), 
                                  x = x_ges$x)
  list_df_10[[i]] <- generate_data(seeds[i], n = n, p = 0.1, b = exp(0.5), 
                                   x = x_ges$x)
  list_df_25[[i]] <- generate_data(seeds[i], n = n, p = 0.25, b = exp(0.5), 
                                   x = x_ges$x)
  list_df_50[[i]] <- generate_data(seeds[i], n = n, p = 0.5, b = exp(0.5), 
                                   x = x_ges$x)
}

## Evaluation data frame:
x_eval <- seq(from = 0, to = 2, length = 1000)
Gx_eval <- expand.grid(xvar = x_eval, G = factor(c("G1", "G2"), levels = c("G1", "G2")))
## Simulate new observations 
set.seed(937) 
sT_G1_test <- sT_G2_test <- rep(0, 1000)
for(i in 1:1000){
  sT_G1_test[i] <- rweibull(1, scale = exp(0.5), shape = (2 + x_eval[i]^2))
  sT_G2_test[i] <- rweibull(1, scale = exp(0.5), shape = (2.5 + x_eval[i]^2))
}

df_eval <- data.frame(sT = c(sT_G1_test, sT_G2_test), G = Gx_eval$G, 
                      xvar = Gx_eval$xvar) 
ngrid_eval_G1 <- sort(unique(df_eval[df_eval$G == "G1",]$sT))
ngrid_eval_G2 <- sort(unique(df_eval[df_eval$G == "G2",]$sT))




##************ Alternative strategies *********************************


## ********************** Cox

list_probs_cox5_G1 <- vector(mode = "list", length = 100)
list_probs_cox10_G1 <- vector(mode = "list", length = 100)
list_probs_cox25_G1 <- vector(mode = "list", length = 100)
list_probs_cox50_G1 <- vector(mode = "list", length = 100)
list_probs_cox5_G2 <- vector(mode = "list", length = 100)
list_probs_cox10_G2 <- vector(mode = "list", length = 100)
list_probs_cox25_G2 <- vector(mode = "list", length = 100)
list_probs_cox50_G2 <- vector(mode = "list", length = 100)

## Estimate Cox models and predict survival probabilites for observations in 
## df_eval over ngrid_eval_G1/ngrid_eval_G2:
for(i in 1:100){

  cox_mod5 <- coxph(Surv(sT, event) ~ G + xvar, data = list_df_5[[i]])
  surv5 <- survfit(cox_mod5, newdata = df_eval)
  surv_probs5 <- matrix(surv5$surv, ncol = 600, byrow = TRUE)
  probs_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  probs_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(k in 1:1000){
    probs_mat_G1[k,] <- stepfun(surv5$time, c(1, surv_probs5[k,]))(ngrid_eval_G1)
    probs_mat_G2[k,] <- stepfun(surv5$time, c(1, surv_probs5[k+1000,]))(ngrid_eval_G2)
  }  
  list_probs_cox5_G1[[i]] <- probs_mat_G1
  list_probs_cox5_G2[[i]] <- probs_mat_G2


  cox_mod10 <- coxph(Surv(sT, event) ~ G + xvar, data = list_df_10[[i]])
  surv10 <- survfit(cox_mod10, newdata = df_eval)
  surv_probs10 <- matrix(surv10$surv, ncol = 600, byrow = TRUE)
  probs_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  probs_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(k in 1:1000){
    probs_mat_G1[k,] <- stepfun(surv10$time, c(1, surv_probs10[k,]))(ngrid_eval_G1)
    probs_mat_G2[k,] <- stepfun(surv10$time, c(1, surv_probs10[k+1000,]))(ngrid_eval_G2)
  }  
  list_probs_cox10_G1[[i]] <- probs_mat_G1
  list_probs_cox10_G2[[i]] <- probs_mat_G2

  cox_mod25 <- coxph(Surv(sT, event) ~ G + xvar, data = list_df_25[[i]])
  surv25 <- survfit(cox_mod25, newdata = df_eval)
  surv_probs25 <- matrix(surv25$surv, ncol = 600, byrow = TRUE)
  probs_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  probs_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(k in 1:1000){
    probs_mat_G1[k,] <- stepfun(surv25$time, c(1, surv_probs25[k,]))(ngrid_eval_G1)
    probs_mat_G2[k,] <- stepfun(surv25$time, c(1, surv_probs25[k+1000,]))(ngrid_eval_G2)
  }  
  list_probs_cox25_G1[[i]] <- probs_mat_G1
  list_probs_cox25_G2[[i]] <- probs_mat_G2

  cox_mod50 <- coxph(Surv(sT, event) ~ G + xvar, data = list_df_50[[i]])
  surv50 <- survfit(cox_mod50, newdata = df_eval)
  surv_probs50 <- matrix(surv50$surv, ncol = 600, byrow = TRUE)
  probs_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  probs_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(k in 1:1000){
    probs_mat_G1[k,] <- stepfun(surv50$time, c(1, surv_probs50[k,]))(ngrid_eval_G1)
    probs_mat_G2[k,] <- stepfun(surv50$time, c(1, surv_probs50[k+1000,]))(ngrid_eval_G2)
  }  
  list_probs_cox50_G1[[i]] <- probs_mat_G1
  list_probs_cox50_G2[[i]] <- probs_mat_G2
}   




## ******************** Stratified Cox

list_probs_stratcox5_G1 <- vector(mode = "list", length = 100)
list_probs_stratcox10_G1 <- vector(mode = "list", length = 100)
list_probs_stratcox25_G1 <- vector(mode = "list", length = 100)
list_probs_stratcox50_G1 <- vector(mode = "list", length = 100)
list_probs_stratcox5_G2 <- vector(mode = "list", length = 100)
list_probs_stratcox10_G2 <- vector(mode = "list", length = 100)
list_probs_stratcox25_G2 <- vector(mode = "list", length = 100)
list_probs_stratcox50_G2 <- vector(mode = "list", length = 100)

## Estimate stratified Cox models and predict survival probabilites for 
## observations in df_eval over ngrid_eval_G1/ngrid_eval_G2:
for(i in 1:100){

  cox_mod5.2 <- coxph(Surv(sT, event) ~ xvar + strata(G), data = list_df_5[[i]])
  surv_5.2 <- survfit(cox_mod5.2, newdata = df_eval)
  surv_mat <- matrix(surv_5.2$surv, ncol = 300, byrow = TRUE)
  times_mat <- matrix(surv_5.2$time, ncol = 300, byrow = TRUE)
  pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(j in 1:1000){
    pred_mat_G1[j,] <- stepfun(times_mat[j,], c(1, surv_mat[j,]))(ngrid_eval_G1)
    pred_mat_G2[j,] <- stepfun(times_mat[j+1000,], c(1, surv_mat[j+1000,]))(ngrid_eval_G2)
  }
  list_probs_stratcox5_G1[[i]] <- pred_mat_G1
  list_probs_stratcox5_G2[[i]] <- pred_mat_G2


  cox_mod10.2 <- coxph(Surv(sT, event) ~ xvar + strata(G), data = list_df_10[[i]])
  surv_10.2 <- survfit(cox_mod10.2, newdata = df_eval)
  surv_mat <- matrix(surv_10.2$surv, ncol = 300, byrow = TRUE)
  times_mat <- matrix(surv_10.2$time, ncol = 300, byrow = TRUE)
  pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(j in 1:1000){
    pred_mat_G1[j,] <- stepfun(times_mat[j,], c(1, surv_mat[j,]))(ngrid_eval_G1)
    pred_mat_G2[j,] <- stepfun(times_mat[j+1000,], c(1, surv_mat[j+1000,]))(ngrid_eval_G2)
  }
  list_probs_stratcox10_G1[[i]] <- pred_mat_G1
  list_probs_stratcox10_G2[[i]] <- pred_mat_G2

  cox_mod25.2 <- coxph(Surv(sT, event) ~ xvar + strata(G), data = list_df_25[[i]])
  surv_25.2 <- survfit(cox_mod25.2, newdata = df_eval)
  surv_mat <- matrix(surv_25.2$surv, ncol = 300, byrow = TRUE)
  times_mat <- matrix(surv_25.2$time, ncol = 300, byrow = TRUE)
  pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(j in 1:1000){
    pred_mat_G1[j,] <- stepfun(times_mat[j,], c(1, surv_mat[j,]))(ngrid_eval_G1)
    pred_mat_G2[j,] <- stepfun(times_mat[j+1000,], c(1, surv_mat[j+1000,]))(ngrid_eval_G2)
  }
  list_probs_stratcox25_G1[[i]] <- pred_mat_G1
  list_probs_stratcox25_G2[[i]] <- pred_mat_G2

  cox_mod50.2 <- coxph(Surv(sT, event) ~ xvar + strata(G), data = list_df_50[[i]])
  surv_50.2 <- survfit(cox_mod50.2, newdata = df_eval)
  surv_mat <- matrix(surv_50.2$surv, ncol = 300, byrow = TRUE)
  times_mat <- matrix(surv_50.2$time, ncol = 300, byrow = TRUE)
  pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
  pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
  for(j in 1:1000){
    pred_mat_G1[j,] <- stepfun(times_mat[j,], c(1, surv_mat[j,]))(ngrid_eval_G1)
    pred_mat_G2[j,] <- stepfun(times_mat[j+1000,], c(1, surv_mat[j+1000,]))(ngrid_eval_G2)
  }
  list_probs_stratcox50_G1[[i]] <- pred_mat_G1
  list_probs_stratcox50_G2[[i]] <- pred_mat_G2
}   



## ******************** Cforest ***************************************

list_probs_cf5_G1 <- list_probs_cf10_G1 <- list_probs_cf25_G1 <- list_probs_cf50_G1 <- vector(mode = "list", length = 100)
list_probs_cf5_G2 <- list_probs_cf10_G2 <- list_probs_cf25_G2 <- list_probs_cf50_G2 <- vector(mode = "list", length = 100)

## Estimate conditional random forests and predict survival probabilities for 
## observations in df_eval over ngrid_eval_G1/ngrid_eval_G2:
for(i in 1:100){
   pred <- predict(cforest(Surv(sT, event) ~ G + xvar, data = list_df_5[[i]], 
                           controls = cforest_control(mtry = 2)),
                   newdata = df_eval, type = "prob")
   pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
   pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
   for(j in 1:1000){
      pred_mat_G1[j,] <- stepfun(pred[[j]]$time, c(1, pred[[j]]$surv))(ngrid_eval_G1)
      pred_mat_G2[j,] <- stepfun(pred[[j+1000]]$time, c(1, pred[[j+1000]]$surv))(ngrid_eval_G2)
   }
   list_probs_cf5_G1[[i]] <- pred_mat_G1 
   list_probs_cf5_G2[[i]] <- pred_mat_G2 

   pred <- predict(cforest(Surv(sT, event) ~ G + xvar, data = list_df_10[[i]], 
                           controls = cforest_control(mtry = 2)),
                   newdata = df_eval, type = "prob")
   pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
   pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
   for(j in 1:1000){
      pred_mat_G1[j,] <- stepfun(pred[[j]]$time, c(1, pred[[j]]$surv))(ngrid_eval_G1)
      pred_mat_G2[j,] <- stepfun(pred[[j+1000]]$time, c(1, pred[[j+1000]]$surv))(ngrid_eval_G2)
   }
   list_probs_cf10_G1[[i]] <- pred_mat_G1 
   list_probs_cf10_G2[[i]] <- pred_mat_G2 

   pred <- predict(cforest(Surv(sT, event) ~ G + xvar, data = list_df_25[[i]], 
                           controls = cforest_control(mtry = 2)),
                   newdata = df_eval, type = "prob")
   pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
   pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
   for(j in 1:1000){
      pred_mat_G1[j,] <- stepfun(pred[[j]]$time, c(1, pred[[j]]$surv))(ngrid_eval_G1)
      pred_mat_G2[j,] <- stepfun(pred[[j+1000]]$time, c(1, pred[[j+1000]]$surv))(ngrid_eval_G2)
   }
   list_probs_cf25_G1[[i]] <- pred_mat_G1 
   list_probs_cf25_G2[[i]] <- pred_mat_G2 

   pred <- predict(cforest(Surv(sT, event) ~ G + xvar, data = list_df_50[[i]], 
                           controls = cforest_control(mtry = 2)),
                   newdata = df_eval, type = "prob")
   pred_mat_G1 <- matrix(0, nrow = 1000, ncol = 1000)
   pred_mat_G2 <- matrix(0, nrow = 1000, ncol = 1000)
   for(j in 1:1000){
      pred_mat_G1[j,] <- stepfun(pred[[j]]$time, c(1, pred[[j]]$surv))(ngrid_eval_G1)
      pred_mat_G2[j,] <- stepfun(pred[[j+1000]]$time, c(1, pred[[j+1000]]$surv))(ngrid_eval_G2)
   }
   list_probs_cf50_G1[[i]] <- pred_mat_G1 
   list_probs_cf50_G2[[i]] <- pred_mat_G2 
}



##************** Conditional transformation models ******************

## ngrid for estimation consists of all unique event and censoring times:
list_ngrid5 <- list_ngrid10 <- list_ngrid25 <- list_ngrid50 <- 
vector(mode = "list", length = 100)
for(i in 1:100){
  list_ngrid5[[i]] <- sort(unique(list_df_5[[i]]$sT))
  list_ngrid10[[i]] <- sort(unique(list_df_10[[i]]$sT))
  list_ngrid25[[i]] <- sort(unique(list_df_25[[i]]$sT))
  list_ngrid50[[i]] <- sort(unique(list_df_50[[i]]$sT))
}

## Estimate Kaplan-Meier estimates of the censoring distribution
list_G5 <- lapply(list_df_5, FUN = function(x){prodlim(Surv(sT, event) ~ 1, 
                                               data = x, reverse = TRUE)})
list_G10 <- lapply(list_df_10, FUN = function(x){prodlim(Surv(sT, event) ~ 1, 
                                                 data = x, reverse = TRUE)})
list_G25 <- lapply(list_df_25, FUN = function(x){prodlim(Surv(sT, event) ~ 1, 
                                                         data = x, reverse = TRUE)})
list_G50 <- lapply(list_df_50, FUN = function(x){prodlim(Surv(sT, event) ~ 1, 
                                                         data = x, reverse = TRUE)})

## Calculate inverse probability of censoring weights
list_ipcw5 <- list_ipcw10 <- list_ipcw25 <- list_ipcw50 <- vector(mode = "list", 
                                                                  length =100)
for(i in 1:100){
  list_ipcw5[[i]] <- inv_weights_graf(ngrid = list_ngrid5[[i]], 
                                      sT = list_df_5[[i]]$sT,
                                      event = list_df_5[[i]]$event, 
                                      G = list_G5[[i]])
  list_ipcw10[[i]] <- inv_weights_graf(ngrid = list_ngrid10[[i]], 
                                       sT = list_df_10[[i]]$sT, 
                                       event = list_df_10[[i]]$event, 
                                       G = list_G10[[i]])
  list_ipcw25[[i]] <- inv_weights_graf(ngrid = list_ngrid25[[i]], 
                                       sT = list_df_25[[i]]$sT, 
                                       event = list_df_25[[i]]$event, 
                                       G = list_G25[[i]])
  list_ipcw50[[i]] <- inv_weights_graf(ngrid = list_ngrid50[[i]], 
                                       sT = list_df_50[[i]]$sT, 
                                       event = list_df_50[[i]]$event, 
                                       G = list_G50[[i]])
}
save(list_ipcw5, file = "list_ipcw5_sim4.Rda")
save(list_ipcw10, file = "list_ipcw10_sim4.Rda")
save(list_ipcw25, file = "list_ipcw25_sim4.Rda")
save(list_ipcw50, file = "list_ipcw50_sim4.Rda")


## ************************* CTM ************************************************

## Setup model formula:
## A separate smooth function over time is estimated for each treatment group;
## an interaction surface of the covariate x and time is estimated:
## constraint = "increasing" guarantees monotone increasing function estimates
fm <- "bbs(sT, df = 4, constraint = \"increasing\") ~ bols(G, df = 2) + 
                                                      bols(xvar, df = 2)"
fm <- as.formula(fm)

## Estimate CTMs and predict survival probabilities for observations in df_eval
## over ngrid_eval_G1/ngrid_eval_G2:
list_pred_G1_mat_5 <- list_pred_G2_mat_5 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod5 <- mctm(fm, data = list_df_5[[i]], family = Binomial_extreme(), 
                   ngrid = list_ngrid5[[i]], weights = list_ipcw5[[i]], 
                   control = boost_control(nu = 0.5, mstop = 3000,
                   trace = FALSE))
  pred_G1_5 <- predict(ctm_mod5, newdata = df_eval[df_eval$G == "G1",], 
                       y = ngrid_eval_G1, type = "response")
  list_pred_G1_mat_5[[i]] <- matrix(pred_G1_5, 
                                    ncol = length(ngrid_eval_G1))
  pred_G2_5 <- predict(ctm_mod5, newdata = df_eval[df_eval$G == "G2",], 
                       y = ngrid_eval_G2, type = "response")
  list_pred_G2_mat_5[[i]] <- matrix(pred_G2_5, 
                                    ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_mat_5, file = "list_pred_G1_mat_5_sim4.Rda")
save(list_pred_G2_mat_5, file = "list_pred_G2_mat_5_sim4.Rda")


list_pred_G1_mat_10 <- list_pred_G2_mat_10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod10 <- mctm(fm, data = list_df_10[[i]], family = Binomial_extreme(), 
                    ngrid = list_ngrid10[[i]], weights = list_ipcw10[[i]], 
                    control = boost_control(nu = 0.5, mstop = 3000,
                    trace = FALSE))
  pred_G1_10 <- predict(ctm_mod10, newdata = df_eval[df_eval$G == "G1",], 
                        y = ngrid_eval_G1, type = "response")
  list_pred_G1_mat_10[[i]] <- matrix(pred_G1_10, 
                                     ncol = length(ngrid_eval_G1))
  pred_G2_10 <- predict(ctm_mod10, newdata = df_eval[df_eval$G == "G2",], 
                        y = ngrid_eval_G2, type = "response")
  list_pred_G2_mat_10[[i]] <- matrix(pred_G2_10, 
                                     ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_mat_10, file = "list_pred_G1_mat_10_sim4.Rda")
save(list_pred_G2_mat_10, file = "list_pred_G2_mat_10_sim4.Rda")

list_pred_G1_mat_25 <- list_pred_G2_mat_25 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod25 <- mctm(fm, data = list_df_25[[i]], family = Binomial_extreme(), 
                    ngrid = list_ngrid25[[i]], weights = list_ipcw25[[i]], 
                    control = boost_control(nu = 0.5, mstop = 3000,
                    trace = FALSE))
  pred_G1_25 <- predict(ctm_mod25, newdata = df_eval[df_eval$G == "G1",], 
                        y = ngrid_eval_G1, type = "response")
  list_pred_G1_mat_25[[i]] <- matrix(pred_G1_25, 
                                     ncol = length(ngrid_eval_G1))
  pred_G2_25 <- predict(ctm_mod25, newdata = df_eval[df_eval$G == "G2",], 
                        y = ngrid_eval_G2, type = "response")
  list_pred_G2_mat_25[[i]] <- matrix(pred_G2_25, 
                                     ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_mat_25, file = "list_pred_G1_mat_25_sim4.Rda")
save(list_pred_G2_mat_25, file = "list_pred_G2_mat_25_sim4.Rda")

list_pred_G1_mat_50 <- list_pred_G2_mat_50 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod50 <- mctm(fm, data = list_df_50[[i]], family = Binomial_extreme(), 
                    ngrid = list_ngrid50[[i]], weights = list_ipcw50[[i]], 
                    control = boost_control(nu = 0.5, mstop = 3000,
                    trace = FALSE))
  pred_G1_50 <- predict(ctm_mod50, newdata = df_eval[df_eval$G == "G1",], 
                        y = ngrid_eval_G1, type = "response")
  list_pred_G1_mat_50[[i]] <- matrix(pred_G1_50, 
                                     ncol = length(ngrid_eval_G1))
  pred_G2_50 <- predict(ctm_mod50, newdata = df_eval[df_eval$G == "G2",], 
                        y = ngrid_eval_G2, type = "response")
  list_pred_G2_mat_50[[i]] <- matrix(pred_G2_50, 
                                     ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_mat_50, file = "list_pred_G1_mat_50_sim4.Rda")
save(list_pred_G2_mat_50, file = "list_pred_G2_mat_50_sim4.Rda")


## ************************************* CLTM ******************************** 

source("mctm_cltm.R")
fm2 <- "bbs(sT, df = 4, constraint = \"increasing\") ~ NULL"
fm2 <- as.formula(fm2)

## Estimate CLTMs and predict survival probabilities for observations in 
## df_eval over ngrid_eval_G1/ngrid_eval_G2:
list_pred_G1_CLTM_5 <- list_pred_G2_CLTM_5 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod5 <- mctm_cltm(fm2, yname = "sT", 
                        constant = "bols(xvar, df = 4) + bols(G, df = 4)", 
                        varying = "bbs(sT, df = 4, constraint = \"increasing\")", 
                        no.variance = TRUE, data = list_df_5[[i]], 
                        family = Binomial_extreme(), ngrid = list_ngrid5[[i]], 
                        weights = list_ipcw5[[i]], control = boost_control(nu = 0.5, 
                        mstop = 3000, trace = FALSE))   
  pred_G1_5 <- predict(ctm_mod5, newdata = df_eval[df_eval$G == "G1",], 
                       y = ngrid_eval_G1, type = "response")
  list_pred_G1_CLTM_5[[i]] <- matrix(pred_G1_5, 
                                     ncol = length(ngrid_eval_G1))
  pred_G2_5 <- predict(ctm_mod5, newdata = df_eval[df_eval$G == "G2",], 
                       y = ngrid_eval_G2, type = "response")
  list_pred_G2_CLTM_5[[i]] <- matrix(pred_G2_5, 
                                     ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_CLTM_5, file = "list_pred_G1_CLTM_5_sim4.Rda")
save(list_pred_G2_CLTM_5, file = "list_pred_G2_CLTM_5_sim4.Rda")


list_pred_G1_CLTM_10 <- list_pred_G2_CLTM_10 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod10 <- mctm_cltm(fm2, yname = "sT", 
                         constant = "bols(xvar, df = 4) + bols(G, df = 4)", 
                         varying = "bbs(sT, df = 4, constraint = \"increasing\")", 
                         no.variance = TRUE, data = list_df_10[[i]], 
                         family = Binomial_extreme(), ngrid = list_ngrid10[[i]], 
                         weights = list_ipcw10[[i]], control = boost_control(nu = 0.5, 
                         mstop = 3000, trace = FALSE))   
  pred_G1_10 <- predict(ctm_mod10, newdata = df_eval[df_eval$G == "G1",], 
                        y = ngrid_eval_G1, type = "response")
  list_pred_G1_CLTM_10[[i]] <- matrix(pred_G1_10, 
                                      ncol = length(ngrid_eval_G1))
  pred_G2_10 <- predict(ctm_mod10, newdata = df_eval[df_eval$G == "G2",], 
                        y = ngrid_eval_G2, type = "response")
  list_pred_G2_CLTM_10[[i]] <- matrix(pred_G2_10, 
                                      ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_CLTM_10, file = "list_pred_G1_CLTM_10_sim4.Rda")
save(list_pred_G2_CLTM_10, file = "list_pred_G2_CLTM_10_sim4.Rda")

list_pred_G1_CLTM_25 <- list_pred_G2_CLTM_25 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod25 <- mctm_cltm(fm2, yname = "sT", 
                         constant = "bols(xvar, df = 4) + bols(G, df = 4)", 
                         varying = "bbs(sT, df = 4, constraint = \"increasing\")", 
                         no.variance = TRUE, data = list_df_25[[i]], 
                         family = Binomial_extreme(), ngrid = list_ngrid25[[i]], 
                         weights = list_ipcw25[[i]], control = boost_control(nu = 0.5, 
                         mstop = 3000, trace = FALSE))  
  pred_G1_25 <- predict(ctm_mod25, newdata = df_eval[df_eval$G == "G1",], 
                        y = ngrid_eval_G1, type = "response")
  list_pred_G1_CLTM_25[[i]] <- matrix(pred_G1_25, 
                                      ncol = length(ngrid_eval_G1))
  pred_G2_25 <- predict(ctm_mod25, newdata = df_eval[df_eval$G == "G2",], 
                        y = ngrid_eval_G2, type = "response")
  list_pred_G2_CLTM_25[[i]] <- matrix(pred_G2_25, 
                                      ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_CLTM_25, file = "list_pred_G1_CLTM_25_sim4.Rda")
save(list_pred_G2_CLTM_25, file = "list_pred_G2_CLTM_25_sim4.Rda")

list_pred_G1_CLTM_50 <- list_pred_G2_CLTM_50 <- vector(mode = "list", length = 100)
for(i in 1:100){
  ctm_mod50 <- mctm_cltm(fm2, yname = "sT", 
                         constant = "bols(xvar, df = 4) + bols(G, df = 4)", 
                         varying = "bbs(sT, df = 4, constraint = \"increasing\")", 
                         no.variance = TRUE, data = list_df_50[[i]], 
                         family = Binomial_extreme(), ngrid = list_ngrid50[[i]], 
                         weights = list_ipcw50[[i]], control = boost_control(nu = 0.5, 
                         mstop = 3000, trace = FALSE))  
  pred_G1_50 <- predict(ctm_mod50, newdata = df_eval[df_eval$G == "G1",], 
                        y = ngrid_eval_G1, type = "response")
  list_pred_G1_CLTM_50[[i]] <- matrix(pred_G1_50, 
                                      ncol = length(ngrid_eval_G1))
  pred_G2_50 <- predict(ctm_mod50, newdata = df_eval[df_eval$G == "G2",], 
                        y = ngrid_eval_G2, type = "response")
  list_pred_G2_CLTM_50[[i]] <- matrix(pred_G2_50, 
                                      ncol = length(ngrid_eval_G2))
}
save(list_pred_G1_CLTM_50, file = "list_pred_G1_CLTM_50_sim4.Rda")
save(list_pred_G2_CLTM_50, file = "list_pred_G2_CLTM_50_sim4.Rda")





##********************* Evaluation **********************************

## Calculate MAD
## N = 2000
## ngrid_G1 / ngrid_G2 = ngrid_eval_G1 / ngrid_eval_G2
## surv_probs_G1/surv_probs_G2: predicted survival probabilites for 
##                              the treatment groups.

MAD <- function(surv_probs_G1, surv_probs_G2, ngrid_G1, ngrid_G2, N, df_eval,
                model_name){
  sum_G1 <- 0 
  for(l in 1:(N/2)){
    c <- 2 + df_eval$x[l]^2
    sum_G1 <- sum_G1 + sum(abs(surv_probs_G1[l,] - Strue(ngrid_G1, b = exp(0.5), 
                                                         c = c)))/length(ngrid_G1)
  }
  MAD_G1 <- sum_G1/(N/2) ## MAD per observation AND grid point

  sum_G2 <- 0
  for(l in 1:(N/2)){
    c <- 2.5 + df_eval$x[l]^2
    sum_G2 <- sum_G2 + sum(abs(surv_probs_G2[l,] - Strue(ngrid_G2, b = exp(0.5), 
                                                         c = c)))/length(ngrid_G2)
  }
  MAD_G2 <- sum_G2/(N/2) ## MAD per observation AND grid point
 
  MAD <- matrix(c(MAD_G1, MAD_G2), nrow = 1, ncol = 2)
  rownames(MAD) <- model_name
  colnames(MAD) <- c("G1", "G2")

  return(MAD)
}


## Calculate log score

## sT_G1_test = df_eval[df_eval$G == "G1",]$sT; length(sT_G1_test) = 1000.
## sT_G2_test = df_eval[df_eval$G == "G2",]$sT; length(sT_G2_test) = 1000.
## = new survival times; 
## ngrid_G1 = ngrid_eval_G1; length(ngrid) = 1000.
## ngrid_G2 = ngrid_eval_G2; length(ngrid) = 1000.
## surv_probs_G1/surv_probs_G2: predicted survival probabilites for 
##                              the treatment groups.

log_score_eval <- function(sT_G1_test, sT_G2_test, surv_probs_G1, surv_probs_G2,
                           ngrid_G1, ngrid_G2, model_name){

  pdf_probs_G1 <- 1 - surv_probs_G1
  pdf_probs_G2 <- 1 - surv_probs_G2

  ls_per_obs_G1 <- 0
  for(j in 1:length(sT_G1_test)){
   ls_per_obs_G1 <- ls_per_obs_G1 + log_score(sT_G1_test[j], ngrid_G1, pdf_probs_G1[j,])
  }
  ls_G1 <- ls_per_obs_G1/length(sT_G1_test)

  ls_per_obs_G2 <- 0
  for(j in 1:length(sT_G2_test)){
   ls_per_obs_G2 <- ls_per_obs_G2 + log_score(sT_G2_test[j], ngrid_G2, pdf_probs_G2[j,])
  }
  ls_G2 <- ls_per_obs_G2/length(sT_G2_test)

  ls <- matrix(c(ls_G1, ls_G2), nrow = 1, ncol = 2)
  rownames(ls) <- model_name
  colnames(ls) <- c("G1", "G2")

  return(ls)
}



## ********************************** MADs *****************************************


## ****************************** Cox

list_MADs_cox5 <- list_MADs_cox10 <- list_MADs_cox25 <- list_MADs_cox50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){
list_MADs_cox5[[i]] <- MAD(list_probs_cox5_G1[[i]], 
                           list_probs_cox5_G2[[i]], 
                           ngrid_G1 = ngrid_eval_G1,
                           ngrid_G2 = ngrid_eval_G2, N = 2000,
                           df_eval, model_name = "Cox")
list_MADs_cox10[[i]] <- MAD(list_probs_cox10_G1[[i]], 
                            list_probs_cox10_G2[[i]], 
                            ngrid_G1 = ngrid_eval_G1, 
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "Cox")
list_MADs_cox25[[i]] <- MAD(list_probs_cox25_G1[[i]], 
                            list_probs_cox25_G2[[i]], 
                            ngrid_G1 = ngrid_eval_G1, 
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "Cox")
list_MADs_cox50[[i]] <- MAD(list_probs_cox50_G1[[i]], 
                            list_probs_cox50_G2[[i]], 
                            ngrid_G1 = ngrid_eval_G1,
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "Cox")
}

MAD_cox5 <- matrix(unlist(list_MADs_cox5), nrow = 2)
rownames(MAD_cox5) <- c("G1", "G2")
mean_MAD_cox5 <- rowMeans(MAD_cox5)
MAD_cox10 <- matrix(unlist(list_MADs_cox10), nrow = 2)
rownames(MAD_cox10) <- c("G1", "G2")
mean_MAD_cox10 <- rowMeans(MAD_cox10)
MAD_cox25 <- matrix(unlist(list_MADs_cox25), nrow = 2)
rownames(MAD_cox25) <- c("G1", "G2")
mean_MAD_cox25 <- rowMeans(MAD_cox25)
MAD_cox50 <- matrix(unlist(list_MADs_cox50), nrow = 2)
rownames(MAD_cox50) <- c("G1", "G2")
mean_MAD_cox50 <- rowMeans(MAD_cox50)
## Save MADs for box plots:
save(MAD_cox5, MAD_cox10, MAD_cox25, MAD_cox50, 
     file = "MADs_cox_sim4_boxplot.Rda")

mean_MAD_cox <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_MAD_cox) <- c("G1", "G2")
colnames(mean_MAD_cox) <- c("5% censoring", "10% censoring", "25% censoring", 
                            "50% censoring") 
mean_MAD_cox[,1] <- mean_MAD_cox5
mean_MAD_cox[,2] <- mean_MAD_cox10
mean_MAD_cox[,3] <- mean_MAD_cox25
mean_MAD_cox[,4] <- mean_MAD_cox50
save(mean_MAD_cox, file = "mean_MAD_cox_sim4.Rda")


## ****************************** Stratified Cox

list_MADs_stratcox5 <- list_MADs_stratcox10 <- list_MADs_stratcox25 <- list_MADs_stratcox50 <- vector(mode = "list", length = 100)

for(i in 1:100){
list_MADs_stratcox5[[i]] <- MAD(list_probs_stratcox5_G1[[i]], 
                                list_probs_stratcox5_G2[[i]], 
                                ngrid_G1 = ngrid_eval_G1,
                                ngrid_G2 = ngrid_eval_G2, N = 2000,
                                df_eval, model_name = "StratCox")
list_MADs_stratcox10[[i]] <- MAD(list_probs_stratcox10_G1[[i]], 
                                 list_probs_stratcox10_G2[[i]], 
                                 ngrid_G1 = ngrid_eval_G1, 
                                 ngrid_G2 = ngrid_eval_G2, N = 2000,
                                 df_eval, model_name = "StratCox")
list_MADs_stratcox25[[i]] <- MAD(list_probs_stratcox25_G1[[i]], 
                                 list_probs_stratcox25_G2[[i]], 
                                 ngrid_G1 = ngrid_eval_G1, 
                                 ngrid_G2 = ngrid_eval_G2, N = 2000,
                                 df_eval, model_name = "StratCox")
list_MADs_stratcox50[[i]] <- MAD(list_probs_stratcox50_G1[[i]], 
                                 list_probs_stratcox50_G2[[i]], 
                                 ngrid_G1 = ngrid_eval_G1,
                                 ngrid_G2 = ngrid_eval_G2, N = 2000,
                                 df_eval, model_name = "StratCox")
}

MAD_stratcox5 <- matrix(unlist(list_MADs_stratcox5), nrow = 2)
rownames(MAD_stratcox5) <- c("G1", "G2")
mean_MAD_stratcox5 <- rowMeans(MAD_stratcox5)
MAD_stratcox10 <- matrix(unlist(list_MADs_stratcox10), nrow = 2)
rownames(MAD_stratcox10) <- c("G1", "G2")
mean_MAD_stratcox10 <- rowMeans(MAD_stratcox10)
MAD_stratcox25 <- matrix(unlist(list_MADs_stratcox25), nrow = 2)
rownames(MAD_stratcox25) <- c("G1", "G2")
mean_MAD_stratcox25 <- rowMeans(MAD_stratcox25)
MAD_stratcox50 <- matrix(unlist(list_MADs_stratcox50), nrow = 2)
rownames(MAD_stratcox50) <- c("G1", "G2")
mean_MAD_stratcox50 <- rowMeans(MAD_stratcox50)
## Save MADs for box plots:
save(MAD_stratcox5, MAD_stratcox10, MAD_stratcox25, MAD_stratcox50, 
     file = "MADs_stratcox_sim4_boxplot.Rda")

mean_MAD_stratcox <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_MAD_stratcox) <- c("G1", "G2")
colnames(mean_MAD_stratcox) <- c("5% censoring", "10% censoring", "25% censoring", 
                                 "50% censoring") 
mean_MAD_stratcox[,1] <- mean_MAD_stratcox5
mean_MAD_stratcox[,2] <- mean_MAD_stratcox10
mean_MAD_stratcox[,3] <- mean_MAD_stratcox25
mean_MAD_stratcox[,4] <- mean_MAD_stratcox50
save(mean_MAD_stratcox, file = "mean_MAD_stratcox_sim4.Rda")



## ****************************** Cforest

list_MADs_cf5 <- list_MADs_cf10 <- list_MADs_cf25 <- list_MADs_cf50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){
list_MADs_cf5[[i]] <- MAD(list_probs_cf5_G1[[i]], 
                          list_probs_cf5_G2[[i]], 
                          ngrid_G1 = ngrid_eval_G1,
                          ngrid_G2 = ngrid_eval_G2, N = 2000,
                          df_eval, model_name = "Cforest")
list_MADs_cf10[[i]] <- MAD(list_probs_cf10_G1[[i]], 
                           list_probs_cf10_G2[[i]], 
                           ngrid_G1 = ngrid_eval_G1,
                           ngrid_G2 = ngrid_eval_G2, N = 2000,
                           df_eval, model_name = "Cforest")
list_MADs_cf25[[i]] <- MAD(list_probs_cf25_G1[[i]], 
                           list_probs_cf25_G2[[i]], 
                           ngrid_G1 = ngrid_eval_G1,
                           ngrid_G2 = ngrid_eval_G2, N = 2000,
                           df_eval, model_name = "Cforest")
list_MADs_cf50[[i]] <- MAD(list_probs_cf50_G1[[i]], 
                           list_probs_cf50_G2[[i]], 
                           ngrid_G1 = ngrid_eval_G1,
                           ngrid_G2 = ngrid_eval_G2, N = 2000,
                           df_eval, model_name = "Cforest")
}

MAD_cf5 <- matrix(unlist(list_MADs_cf5), nrow = 2)
rownames(MAD_cf5) <- c("G1", "G2")
mean_MAD_cf5 <- rowMeans(MAD_cf5)
MAD_cf10 <- matrix(unlist(list_MADs_cf10), nrow = 2)
rownames(MAD_cf10) <- c("G1", "G2")
mean_MAD_cf10 <- rowMeans(MAD_cf10)
MAD_cf25 <- matrix(unlist(list_MADs_cf25), nrow = 2)
rownames(MAD_cf25) <- c("G1", "G2")
mean_MAD_cf25 <- rowMeans(MAD_cf25)
MAD_cf50 <- matrix(unlist(list_MADs_cf50), nrow = 2)
rownames(MAD_cf50) <- c("G1", "G2")
mean_MAD_cf50 <- rowMeans(MAD_cf50)
## Save MADs for box plots:
save(MAD_cf5, MAD_cf10, MAD_cf25, MAD_cf50, 
     file = "MADs_cf_sim4_boxplot.Rda")

mean_MAD_cf <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_MAD_cf) <- c("G1", "G2")
colnames(mean_MAD_cf) <- c("5% censoring", "10% censoring", "25% censoring", 
                           "50% censoring") 
mean_MAD_cf[,1] <- mean_MAD_cf5
mean_MAD_cf[,2] <- mean_MAD_cf10
mean_MAD_cf[,3] <- mean_MAD_cf25
mean_MAD_cf[,4] <- mean_MAD_cf50
save(mean_MAD_cf, file = "mean_MAD_cf_sim4.Rda")



## ****************************** CTM

list_MADs_CTM5 <- list_MADs_CTM10 <- list_MADs_CTM25 <- list_MADs_CTM50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){
list_MADs_CTM5[[i]] <- MAD(1 - list_pred_G1_mat_5[[i]], 
                           1 - list_pred_G2_mat_5[[i]], 
                           ngrid_G1 = ngrid_eval_G1, 
                           ngrid_G2 = ngrid_eval_G2, N = 2000,
                           df_eval, model_name = "CTM")
list_MADs_CTM10[[i]] <- MAD(1 - list_pred_G1_mat_10[[i]], 
                            1 - list_pred_G2_mat_10[[i]], 
                            ngrid_G1 = ngrid_eval_G1, 
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "CTM")
list_MADs_CTM25[[i]] <- MAD(1 - list_pred_G1_mat_25[[i]], 
                            1 - list_pred_G2_mat_25[[i]], 
                            ngrid_G1 = ngrid_eval_G1, 
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "CTM")
list_MADs_CTM50[[i]] <- MAD(1 - list_pred_G1_mat_50[[i]], 
                            1 - list_pred_G2_mat_50[[i]], 
                            ngrid_G1 = ngrid_eval_G1, 
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "CTM")
}

MAD_CTM5 <- matrix(unlist(list_MADs_CTM5), nrow = 2)
rownames(MAD_CTM5) <- c("G1", "G2")
mean_MAD_CTM5 <- rowMeans(MAD_CTM5)
MAD_CTM10 <- matrix(unlist(list_MADs_CTM10), nrow = 2)
rownames(MAD_CTM10) <- c("G1", "G2")
mean_MAD_CTM10 <- rowMeans(MAD_CTM10)
MAD_CTM25 <- matrix(unlist(list_MADs_CTM25), nrow = 2)
rownames(MAD_CTM25) <- c("G1", "G2")
mean_MAD_CTM25 <- rowMeans(MAD_CTM25)
MAD_CTM50 <- matrix(unlist(list_MADs_CTM50), nrow = 2)
rownames(MAD_CTM50) <- c("G1", "G2")
mean_MAD_CTM50 <- rowMeans(MAD_CTM50)
## Save MADs for box plots:
save(MAD_CTM5, MAD_CTM10, MAD_CTM25, MAD_CTM50, 
     file = "MADs_CTM_sim4_boxplot.Rda")

mean_MAD_CTM <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_MAD_CTM) <- c("G1", "G2")
colnames(mean_MAD_CTM) <- c("5% censoring", "10% censoring", "25% censoring", 
                            "50% censoring") 
mean_MAD_CTM[,1] <- mean_MAD_CTM5
mean_MAD_CTM[,2] <- mean_MAD_CTM10
mean_MAD_CTM[,3] <- mean_MAD_CTM25
mean_MAD_CTM[,4] <- mean_MAD_CTM50
save(mean_MAD_CTM, file = "mean_MAD_CTM_sim4.Rda")



## ****************************** CLTM

list_MADs_CLTM5 <- list_MADs_CLTM10 <- list_MADs_CLTM25 <- list_MADs_CLTM50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){
list_MADs_CLTM5[[i]] <- MAD(1 - list_pred_G1_CLTM_5[[i]], 
                            1 - list_pred_G2_CLTM_5[[i]], 
                            ngrid_G1 = ngrid_eval_G1,
                            ngrid_G2 = ngrid_eval_G2, N = 2000,
                            df_eval, model_name = "CLTM")
list_MADs_CLTM10[[i]] <- MAD(1 - list_pred_G1_CLTM_10[[i]], 
                             1 - list_pred_G2_CLTM_10[[i]], 
                             ngrid_G1 = ngrid_eval_G1,
                             ngrid_G2 = ngrid_eval_G2, N = 2000,
                             df_eval, model_name = "CLTM")
list_MADs_CLTM25[[i]] <- MAD(1 - list_pred_G1_CLTM_25[[i]], 
                             1 - list_pred_G2_CLTM_25[[i]], 
                             ngrid_G1 = ngrid_eval_G1,
                             ngrid_G2 = ngrid_eval_G2, N = 2000,
                             df_eval, model_name = "CLTM")
list_MADs_CLTM50[[i]] <- MAD(1 - list_pred_G1_CLTM_50[[i]], 
                             1 - list_pred_G2_CLTM_50[[i]], 
                             ngrid_G1 = ngrid_eval_G1,
                             ngrid_G2 = ngrid_eval_G2, N = 2000,
                             df_eval, model_name = "CLTM")
}

MAD_CLTM5 <- matrix(unlist(list_MADs_CLTM5), nrow = 2)
rownames(MAD_CLTM5) <- c("G1", "G2")
mean_MAD_CLTM5 <- rowMeans(MAD_CLTM5)
MAD_CLTM10 <- matrix(unlist(list_MADs_CLTM10), nrow = 2)
rownames(MAD_CLTM10) <- c("G1", "G2")
mean_MAD_CLTM10 <- rowMeans(MAD_CLTM10)
MAD_CLTM25 <- matrix(unlist(list_MADs_CLTM25), nrow = 2)
rownames(MAD_CLTM25) <- c("G1", "G2")
mean_MAD_CLTM25 <- rowMeans(MAD_CLTM25)
MAD_CLTM50 <- matrix(unlist(list_MADs_CLTM50), nrow = 2)
rownames(MAD_CLTM50) <- c("G1", "G2")
mean_MAD_CLTM50 <- rowMeans(MAD_CLTM50)
## Save MADs for box plots:
save(MAD_CLTM5, MAD_CLTM10, MAD_CLTM25, MAD_CLTM50, 
     file = "MADs_CLTM_sim4_boxplot.Rda")

mean_MAD_CLTM <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_MAD_CLTM) <- c("G1", "G2")
colnames(mean_MAD_CLTM) <- c("5% censoring", "10% censoring", "25% censoring", 
                             "50% censoring") 
mean_MAD_CLTM[,1] <- mean_MAD_CLTM5
mean_MAD_CLTM[,2] <- mean_MAD_CLTM10
mean_MAD_CLTM[,3] <- mean_MAD_CLTM25
mean_MAD_CLTM[,4] <- mean_MAD_CLTM50
save(mean_MAD_CLTM, file = "mean_MAD_CLTM_sim4.Rda")




## ********************************* Log Scores ************************************


## ************************* Cox
list_ls_cox5 <- list_ls_cox10 <- list_ls_cox25 <- list_ls_cox50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){

  list_ls_cox5[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                      df_eval$sT[1001:2000],
                                      list_probs_cox5_G1[[i]], 
                                      list_probs_cox5_G2[[i]],
                                      ngrid_G1 = ngrid_eval_G1,
                                      ngrid_G2 = ngrid_eval_G2, "Cox")
  list_ls_cox10[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       list_probs_cox10_G1[[i]], 
                                       list_probs_cox10_G2[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "Cox")
  list_ls_cox25[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       list_probs_cox25_G1[[i]], 
                                       list_probs_cox25_G2[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "Cox")
  list_ls_cox50[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       list_probs_cox50_G1[[i]], 
                                       list_probs_cox50_G2[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "Cox")
}
ls_cox5 <- matrix(unlist(list_ls_cox5), nrow = 2)
rownames(ls_cox5) <- c("G1", "G2")
mean_ls_cox5 <- rowMeans(ls_cox5)
ls_cox10 <- matrix(unlist(list_ls_cox10), nrow = 2)
rownames(ls_cox10) <- c("G1", "G2")
mean_ls_cox10 <- rowMeans(ls_cox10)
ls_cox25 <- matrix(unlist(list_ls_cox25), nrow = 2)
rownames(ls_cox25) <- c("G1", "G2")
mean_ls_cox25 <- rowMeans(ls_cox25)
ls_cox50 <- matrix(unlist(list_ls_cox50), nrow = 2)
rownames(ls_cox50) <- c("G1", "G2")
mean_ls_cox50 <- rowMeans(ls_cox50)
## save log scores for boxplots
save(ls_cox5, ls_cox10, ls_cox25, ls_cox50, 
     file = "ls_cox_sim4_boxplot.Rda") 

mean_ls_cox <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_ls_cox) <- c("G1", "G2")
colnames(mean_ls_cox) <- c("5% censoring", "10% censoring", "25% censoring", 
                      "50% censoring") 
mean_ls_cox[,1] <- mean_ls_cox5
mean_ls_cox[,2] <- mean_ls_cox10
mean_ls_cox[,3] <- mean_ls_cox25
mean_ls_cox[,4] <- mean_ls_cox50

save(mean_ls_cox, file = "mean_ls_cox_sim4.Rda")



## ************************* Stratified Cox
list_ls_stratcox5 <- list_ls_stratcox10 <- list_ls_stratcox25 <- list_ls_stratcox50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){

  list_ls_stratcox5[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                           df_eval$sT[1001:2000],
                                           list_probs_stratcox5_G1[[i]], 
                                           list_probs_stratcox5_G2[[i]],
                                           ngrid_G1 = ngrid_eval_G1,
                                           ngrid_G2 = ngrid_eval_G2, 
                                           "StratCox")
  list_ls_stratcox10[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                            df_eval$sT[1001:2000],
                                            list_probs_stratcox10_G1[[i]], 
                                            list_probs_stratcox10_G2[[i]],
                                            ngrid_G1 = ngrid_eval_G1,
                                            ngrid_G2 = ngrid_eval_G2, 
                                            "StratCox")
  list_ls_stratcox25[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                            df_eval$sT[1001:2000],
                                            list_probs_stratcox25_G1[[i]], 
                                            list_probs_stratcox25_G2[[i]],
                                            ngrid_G1 = ngrid_eval_G1,
                                            ngrid_G2 = ngrid_eval_G2, 
                                            "StratCox")
  list_ls_stratcox50[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                            df_eval$sT[1001:2000],
                                            list_probs_stratcox50_G1[[i]], 
                                            list_probs_stratcox50_G2[[i]],
                                            ngrid_G1 = ngrid_eval_G1,
                                            ngrid_G2 = ngrid_eval_G2, 
                                            "StratCox")
}
ls_stratcox5 <- matrix(unlist(list_ls_stratcox5), nrow = 2)
rownames(ls_stratcox5) <- c("G1", "G2")
mean_ls_stratcox5 <- rowMeans(ls_stratcox5)
ls_stratcox10 <- matrix(unlist(list_ls_stratcox10), nrow = 2)
rownames(ls_stratcox10) <- c("G1", "G2")
mean_ls_stratcox10 <- rowMeans(ls_stratcox10)
ls_stratcox25 <- matrix(unlist(list_ls_stratcox25), nrow = 2)
rownames(ls_stratcox25) <- c("G1", "G2")
mean_ls_stratcox25 <- rowMeans(ls_stratcox25)
ls_stratcox50 <- matrix(unlist(list_ls_stratcox50), nrow = 2)
rownames(ls_stratcox50) <- c("G1", "G2")
mean_ls_stratcox50 <- rowMeans(ls_stratcox50)
## save log scores for boxplots
save(ls_stratcox5, ls_stratcox10, ls_stratcox25, ls_stratcox50, 
     file = "ls_stratcox_sim4_boxplot.Rda") 

mean_ls_stratcox <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_ls_stratcox) <- c("G1", "G2")
colnames(mean_ls_stratcox) <- c("5% censoring", "10% censoring", "25% censoring", 
                                "50% censoring") 
mean_ls_stratcox[,1] <- mean_ls_stratcox5
mean_ls_stratcox[,2] <- mean_ls_stratcox10
mean_ls_stratcox[,3] <- mean_ls_stratcox25
mean_ls_stratcox[,4] <- mean_ls_stratcox50

save(mean_ls_stratcox, file = "mean_ls_stratcox_sim4.Rda")



## ************************* Cforest
list_ls_cf5 <- list_ls_cf10 <- list_ls_cf25 <- list_ls_cf50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){

  list_ls_cf5[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                     df_eval$sT[1001:2000],
                                     list_probs_cf5_G1[[i]], 
                                     list_probs_cf5_G2[[i]],
                                     ngrid_G1 = ngrid_eval_G1,
                                     ngrid_G2 = ngrid_eval_G2, "Cforest")
  list_ls_cf10[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                      df_eval$sT[1001:2000],
                                      list_probs_cf10_G1[[i]], 
                                      list_probs_cf10_G2[[i]],
                                      ngrid_G1 = ngrid_eval_G1,
                                      ngrid_G2 = ngrid_eval_G2, "Cforest")
  list_ls_cf25[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                      df_eval$sT[1001:2000],
                                      list_probs_cf25_G1[[i]], 
                                      list_probs_cf25_G2[[i]],
                                      ngrid_G1 = ngrid_eval_G1,
                                      ngrid_G2 = ngrid_eval_G2, "Cforest")
  list_ls_cf50[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                      df_eval$sT[1001:2000],
                                      list_probs_cf50_G1[[i]], 
                                      list_probs_cf50_G2[[i]],
                                      ngrid_G1 = ngrid_eval_G1,
                                      ngrid_G2 = ngrid_eval_G2, "Cforest")
}
ls_cf5 <- matrix(unlist(list_ls_cf5), nrow = 2)
rownames(ls_cf5) <- c("G1", "G2")
mean_ls_cf5 <- rowMeans(ls_cf5)
ls_cf10 <- matrix(unlist(list_ls_cf10), nrow = 2)
rownames(ls_cf10) <- c("G1", "G2")
mean_ls_cf10 <- rowMeans(ls_cf10)
ls_cf25 <- matrix(unlist(list_ls_cf25), nrow = 2)
rownames(ls_cf25) <- c("G1", "G2")
mean_ls_cf25 <- rowMeans(ls_cf25)
ls_cf50 <- matrix(unlist(list_ls_cf50), nrow = 2)
rownames(ls_cf50) <- c("G1", "G2")
mean_ls_cf50 <- rowMeans(ls_cf50)
## save log scores for boxplots
save(ls_cf5, ls_cf10, ls_cf25, ls_cf50, 
     file = "ls_cf_sim4_boxplot.Rda") 

mean_ls_cf <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_ls_cf) <- c("G1", "G2")
colnames(mean_ls_cf) <- c("5% censoring", "10% censoring", "25% censoring", 
                          "50% censoring") 
mean_ls_cf[,1] <- mean_ls_cf5
mean_ls_cf[,2] <- mean_ls_cf10
mean_ls_cf[,3] <- mean_ls_cf25
mean_ls_cf[,4] <- mean_ls_cf50

save(mean_ls_cf, file = "mean_ls_cf_sim4.Rda")


## ************************* CTM
list_ls_CTM5 <- list_ls_CTM10 <- list_ls_CTM25 <- list_ls_CTM50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){

  list_ls_CTM5[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                      df_eval$sT[1001:2000],
                                      1 - list_pred_G1_mat_5[[i]], 
                                      1 - list_pred_G2_mat_5[[i]],
                                      ngrid_G1 = ngrid_eval_G1,
                                      ngrid_G2 = ngrid_eval_G2, "CTM")
  list_ls_CTM10[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       1 - list_pred_G1_mat_10[[i]], 
                                       1 - list_pred_G2_mat_10[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "CTM")
  list_ls_CTM25[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       1 - list_pred_G1_mat_25[[i]], 
                                       1 - list_pred_G2_mat_25[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "CTM")
  list_ls_CTM50[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       1 - list_pred_G1_mat_50[[i]], 
                                       1 - list_pred_G2_mat_50[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "CTM")
}
ls_CTM5 <- matrix(unlist(list_ls_CTM5), nrow = 2)
rownames(ls_CTM5) <- c("G1", "G2")
mean_ls_CTM5 <- rowMeans(ls_CTM5)
ls_CTM10 <- matrix(unlist(list_ls_CTM10), nrow = 2)
rownames(ls_CTM10) <- c("G1", "G2")
mean_ls_CTM10 <- rowMeans(ls_CTM10)
ls_CTM25 <- matrix(unlist(list_ls_CTM25), nrow = 2)
rownames(ls_CTM25) <- c("G1", "G2")
mean_ls_CTM25 <- rowMeans(ls_CTM25)
ls_CTM50 <- matrix(unlist(list_ls_CTM50), nrow = 2)
rownames(ls_CTM50) <- c("G1", "G2")
mean_ls_CTM50 <- rowMeans(ls_CTM50)
## save log scores for boxplots
save(ls_CTM5, ls_CTM10, ls_CTM25, ls_CTM50, 
     file = "ls_CTM_sim4_boxplot.Rda") 

mean_ls_CTM <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_ls_CTM) <- c("G1", "G2")
colnames(mean_ls_CTM) <- c("5% censoring", "10% censoring", "25% censoring", 
                           "50% censoring") 
mean_ls_CTM[,1] <- mean_ls_CTM5
mean_ls_CTM[,2] <- mean_ls_CTM10
mean_ls_CTM[,3] <- mean_ls_CTM25
mean_ls_CTM[,4] <- mean_ls_CTM50

save(mean_ls_CTM, file = "mean_ls_CTM_sim4.Rda")



## ************************* CLTM
list_ls_CLTM5 <- list_ls_CLTM10 <- list_ls_CLTM25 <- list_ls_CLTM50 <- 
vector(mode = "list", length = 100)

for(i in 1:100){

  list_ls_CLTM5[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                       df_eval$sT[1001:2000],
                                       1 - list_pred_G1_CLTM_5[[i]], 
                                       1 - list_pred_G2_CLTM_5[[i]],
                                       ngrid_G1 = ngrid_eval_G1,
                                       ngrid_G2 = ngrid_eval_G2, "CLTM")
  list_ls_CLTM10[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                        df_eval$sT[1001:2000],
                                        1 - list_pred_G1_CLTM_10[[i]], 
                                        1 - list_pred_G2_CLTM_10[[i]],
                                        ngrid_G1 = ngrid_eval_G1,
                                        ngrid_G2 = ngrid_eval_G2, "CLTM")
  list_ls_CLTM25[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                        df_eval$sT[1001:2000],
                                        1 - list_pred_G1_CLTM_25[[i]], 
                                        1 - list_pred_G2_CLTM_25[[i]],
                                        ngrid_G1 = ngrid_eval_G1,
                                        ngrid_G2 = ngrid_eval_G2, "CLTM")
  list_ls_CLTM50[[i]] <- log_score_eval(df_eval$sT[1:1000], 
                                        df_eval$sT[1001:2000],
                                        1 - list_pred_G1_CLTM_50[[i]], 
                                        1 - list_pred_G2_CLTM_50[[i]],
                                        ngrid_G1 = ngrid_eval_G1,
                                        ngrid_G2 = ngrid_eval_G2, "CLTM")
}
ls_CLTM5 <- matrix(unlist(list_ls_CLTM5), nrow = 2)
rownames(ls_CLTM5) <- c("G1", "G2")
mean_ls_CLTM5 <- rowMeans(ls_CLTM5)
ls_CLTM10 <- matrix(unlist(list_ls_CLTM10), nrow = 2)
rownames(ls_CLTM10) <- c("G1", "G2")
mean_ls_CLTM10 <- rowMeans(ls_CLTM10)
ls_CLTM25 <- matrix(unlist(list_ls_CLTM25), nrow = 2)
rownames(ls_CLTM25) <- c("G1", "G2")
mean_ls_CLTM25 <- rowMeans(ls_CLTM25)
ls_CLTM50 <- matrix(unlist(list_ls_CLTM50), nrow = 2)
rownames(ls_CLTM50) <- c("G1", "G2")
mean_ls_CLTM50 <- rowMeans(ls_CLTM50)
## save log scores for boxplots
save(ls_CLTM5, ls_CLTM10, ls_CLTM25, ls_CLTM50, 
     file = "ls_CLTM_sim4_boxplot.Rda") 

mean_ls_CLTM <- matrix(0, nrow = 2, ncol = 4)
rownames(mean_ls_CLTM) <- c("G1", "G2")
colnames(mean_ls_CLTM) <- c("5% censoring", "10% censoring", "25% censoring", 
                            "50% censoring") 
mean_ls_CLTM[,1] <- mean_ls_CLTM5
mean_ls_CLTM[,2] <- mean_ls_CLTM10
mean_ls_CLTM[,3] <- mean_ls_CLTM25
mean_ls_CLTM[,4] <- mean_ls_CLTM50

save(mean_ls_CLTM, file = "mean_ls_CLTM_sim4.Rda")



## ****************************** Tabulars for paper:

## Mean MAD
Mean_MAD <- rbind(mean_MAD_CTM[1,], mean_MAD_CLTM[1,], mean_MAD_cox[1,], 
                  mean_MAD_cf[1,], mean_MAD_stratcox[1,],
                  mean_MAD_CTM[2,], mean_MAD_CLTM[2,], mean_MAD_cox[2,], 
                  mean_MAD_cf[2,], mean_MAD_stratcox[2,])
rownames(Mean_MAD) <- c("CTM G1", "CLTM G1", "Cox G1", "Cforest G1", "StratCox G1", 
                        "CTM G2", "CLTM G2", "Cox G2", "Cforest G2", "StratCox G2")
Mean_MAD <- Mean_MAD * 100
save(Mean_MAD, file = "Mean_MAD_sim4.Rda")

## Mean log score
mean_log_scores <- rbind(mean_ls_CTM[1,], mean_ls_CLTM[1,], mean_ls_cox[1,], 
                         mean_ls_stratcox[1,], mean_ls_cf[1,],
                         mean_ls_CTM[2,], mean_ls_CLTM[2,], mean_ls_cox[2,], 
                         mean_ls_stratcox[2,], mean_ls_cf[2,])
rownames(mean_log_scores) <- c("CTM G1", "CLTM G1", "Cox G1", "StratCox G1", 
                               "Cforest G1",
                               "CTM G2", "CLTM G2", "Cox G2", "StratCox G2", 
                               "Cforest G2")
mean_log_scores <- mean_log_scores * 100
save(mean_log_scores, 
     file = "mean_log_scores_sim4.Rda")
