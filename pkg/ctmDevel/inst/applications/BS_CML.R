##********************************************************************
## Chronic myelogenous leukaemia
##********************************************************************

## The following R Code allows to retrace the analysis of the chronic myelogenous
## leukemia data in Section 4 of the main paper. Since we are not able to deliver
## the original data set, we simulated a comparable data set. The relevant
## explanatory variables and the survival times are sampled randomly and
## independently. The number of observations n = 507 equals the number of
## observations in the original data set. 

## The data set is analysed in terms of five different models:
## 1.) cml_cox: Cox model with linear influences of treatment, risk group, 
##              sex and age.
## 2.) cml_cox_int: The Cox model with interactions includes linear influences
##                  for treatment, risk group, sex and age; all two-times 
##                  interactions between treatment, risk group and sex; 
##                  the three-time interaction between treatment, risk group and
##                  sex.
## 3.) cfr: Conditional random forests with explanatory variables treatment, 
##          risk group, sex and age.
## 4.) cml_ctm: CTM with time-varying influences for treatment, risk group and
##              age.
## 5.) cml_ctm_int: CTM with interactions with time-varying influences for the
##                  treatment - risk group interaction, sex and age.    
## Out-of-sample model performance is measured using the censored log score and
## based on B = 100 bootstrap samples.

## All analysis were carried out in R version 3.0.2. Furthermore, one needs 
## to install the add-on packages ctmDevel (depending on mboostDevel), 
## survival (for the estimation of Cox models), prodlim (for the calculation of
## Kaplan-Meier estimators), and party (for the estimation of conditional
## random forests).

library("ctmDevel")
library("prodlim")
library("survival")
library("party")


## Simulate a comparable data set
set.seed(29)
n <- 507 ## the original data set includes 507 observations
tr <- factor(c("IFN-alpha", "BUS", "HU"), 
             levels = c("IFN-alpha", "BUS", "HU"))
treatment <- sample(tr, n, replace = TRUE)
risk <- factor(c("0", "1", "2"), levels = c("0", "1", "2")) 
riskgroup <- sample(risk, n, replace = TRUE)
se <- factor(c("0", "1"), levels = c("0", "1"))
sex <- sample(se, n, replace = TRUE)
ag <- factor(15:85, levels = 15:85) 
age <- sample(ag, n, replace = TRUE)
age <- as.numeric(age)
ctr <- factor(1:57, levels = 1:57)
center <- sample(ctr, n, replace = TRUE)
## Right-censored event times
T <- round(rexp(n, rate = 0.05), digits = 0) ## event time
C <- round(runif(n, min = 0, max = max(T)), digits = 0) ## censoring time
S <- cbind(T = T, C = C)
sT <- apply(S, 1, min) ## observed event time
cml <- data.frame(treatment = treatment, riskgroup = riskgroup, 
                  sex = sex, age = age, center = center, 
                  time = sT, status = I(T <= C))

cml$riskgroup_treatment <- interaction(cml$treatment, cml$riskgroup)


## ****************** Calculate IPCWs
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



## Function for integrated censored log score
## sT: observed survival time
## ngrid: grid of time points over which survival probabilities were 
##        estimated.
## pis: estimated survival probabilities over ngrid.
## weight: IPCWs over ngrid.
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
  return((-1)*sum_it/length(ngrid))    
}


### Booststrap - data
nboot <- 100  ## number of Bootstrap-Datasets
set.seed(12345)
indx <- vector(mode = "list", length = nboot)
for(i in 1:nboot){
  indx[[i]] <- sample(1:nrow(cml), nrow(cml), replace = TRUE)
}


ngrid <- sort(unique(cml$time))
G <- prodlim(Surv(time, status) ~ 1, data = cml, reverse = TRUE)
ipcw <-   inv_weights_graf(ngrid, sT = cml$time, event = cml$status, G = G)
ipcw_mat <- matrix(ipcw, ncol = length(ngrid))

## CTM
source("Binomial_extreme.R")
## Model formulas
fm <- "bbs(time, df = 2.05, constraint = \"increasing\") ~ 
                     bbs(age, df = 2.05) + 
                     bols(sex, df = 2.05) +
                     bols(riskgroup, df = 2.05) +
                     bols(treatment, df = 2.05)"
fm <- as.formula(fm)
fm_int <- "bbs(time, df = 2.05, constraint = \"increasing\") ~ 
                   bbs(age, df = 2.05) + 
                   bols(sex, df = 2.05) +
                   bols(riskgroup_treatment, df = 2.05)"
fm_int <- as.formula(fm_int)

log_score <- matrix(0, nrow = nboot, ncol = 5)
colnames(log_score) <- c("log score Cox", "log score Cox Int", "log score Cforest",
                         "log score CTM", "log score CTM Int")
rownames(log_score) <- 1:nboot

for(i in 1:nboot){
   print(i)

   data_est <- cml[indx[[i]],]
   data_eval <- cml[-unique(indx[[i]]),]
   indx_eval <- (1:507)[-unique(indx[[i]])]

   ipcw_eval_mat <- matrix(0, nrow = nrow(data_eval), ncol = ncol(ipcw_mat))
   for(l in 1:length(indx_eval)){
     ipcw_eval_mat[l,] <- ipcw_mat[indx_eval[l],]
       
   }

   ## Cox model without interactions
   cml_cox <- coxph(Surv(time, status) ~ treatment + sex + age + riskgroup,  
                    data = data_est)
   surv_cox <- survfit(cml_cox, newdata = data_eval)
   surv_probs <- matrix(surv_cox$surv, ncol = length(unique(data_est$time)), 
                        byrow = TRUE)
   pis_cox <- matrix(0, nrow = nrow(data_eval), ncol = length(ngrid))
   for(k in 1:nrow(data_eval)){
      pis_cox[k,] <- stepfun(surv_cox$time, c(1, surv_probs[k,]))(ngrid)
   }

   pis_cox <- 1 - pis_cox
   ls_temp_cox <- 0
   indx_eval <- as.integer(rownames(data_eval))
   for(l in 1:nrow(data_eval)){
     ls_temp_cox <- ls_temp_cox + log_score_cens(data_eval$time[l], ngrid,
                                                 pis_cox[l,], 
                                                 ipcw_mat[l,])
   }
   log_score[i,1] <- ls_temp_cox/nrow(data_eval)

   ## Cox model with interactions 
   cml_cox_int <- coxph(Surv(time, status) ~ treatment + sex + age + riskgroup +
                    treatment:riskgroup + treatment:sex + riskgroup:sex +
                    treatment:riskgroup:sex,  data = data_est)
   surv_cox <- survfit(cml_cox_int, newdata = data_eval)
   surv_probs <- matrix(surv_cox$surv, ncol = length(unique(data_est$time)), 
                        byrow = TRUE)
   pis_cox <- matrix(0, nrow = nrow(data_eval), ncol = length(ngrid))
   for(k in 1:nrow(data_eval)){
      pis_cox[k,] <- stepfun(surv_cox$time, c(1, surv_probs[k,]))(ngrid)
   }
   pis_cox <- 1 - pis_cox
   ls_temp_cox_int <- 0
   for(l in 1:nrow(data_eval)){
     ls_temp_cox_int <- ls_temp_cox_int + log_score_cens(data_eval$time[l], ngrid,
                                                         pis_cox[l,], 
                                                         ipcw_eval_mat[l,])
   }
   log_score[i,2] <- ls_temp_cox_int/nrow(data_eval)

   ## Cforest
   cfr <- cforest(Surv(time, status) ~ riskgroup + treatment + sex + age, 
                  data = data_est, controls = cforest_control(mtry = 4))

   pred <- predict(cfr, newdata = data_eval, type = "prob")
   pred_cfr <- matrix(0, nrow = nrow(data_eval), ncol = length(ngrid))
   for(k in 1:nrow(data_eval)){
     pred_cfr[k,] <- stepfun(pred[[k]]$time, c(1, pred[[k]]$surv))(ngrid)
   }
   pred_cfr <- 1 - pred_cfr
   ls_temp_cfr <- 0
   for(l in 1:nrow(data_eval)){
     ls_temp_cfr <- ls_temp_cfr + log_score_cens(data_eval$time[l], ngrid,
                                                 pred_cfr[l,], 
                                                 ipcw_eval_mat[l,])
   }
   log_score[i,3] <- ls_temp_cfr/nrow(data_eval)


   ## CTM
   ipcw_est_mat <- matrix(0, nrow = nrow(data_est), ncol = ncol(ipcw_mat))
   for(l in 1:length(indx[[i]])){
     ipcw_est_mat[l,] <- ipcw_mat[indx[[i]][l],]
   }

   options(mboost_useMatrix = FALSE)
   ctm_cml <- mctm(fm, data = data_est, family = Binomial_extreme(),              
                   ngrid = ngrid, weights = as.vector(ipcw_est_mat), 
                   control = boost_control(nu = 0.1, mstop = 1900, 
                                           trace = FALSE))
   ## Optimal mstop = 1900 has been found using multiple resampled data sets.
   ## mstop = 1900 equals the mstop value associated with the smallest
   ## out-of-sample censored log score for new observations.  
   pis <- predict(ctm_cml, newdata = data_eval, y = ngrid, 
                  type = "response")
   pis_mat <- matrix(pis, ncol = length(ngrid))
   ls_temp_CTM <- 0
   for(l in 1:nrow(data_eval)){
     ls_temp_CTM <- ls_temp_CTM + log_score_cens(data_eval$time[l], ngrid, 
                                                 pis_mat[l,], 
                                                 ipcw_eval_mat[l,])
   }
   log_score[i,4] <- ls_temp_CTM/nrow(data_eval)

   ## CTM with interaction
   options(mboost_useMatrix = FALSE)
   ctm_cml_int <- mctm(fm_int, data = data_est, family = Binomial_extreme(),
                       ngrid = ngrid, weights = as.vector(ipcw_est_mat), 
                       control = boost_control(nu = 0.1, mstop = 1900, 
                                               trace = FALSE))
   pis_int <- predict(ctm_cml_int, newdata = data_eval, y = ngrid,
                      type = "response")
   pis_int_mat <- matrix(pis_int, ncol = length(ngrid))
   ls_temp_CTM_int <- 0
   for(l in 1:nrow(data_eval)){
     ls_temp_CTM_int <- ls_temp_CTM_int + log_score_cens(data_eval$time[l], ngrid, 
                                                         pis_int_mat[l,], 
                                                         ipcw_eval_mat[l,])
   }
   log_score[i,5] <- ls_temp_CTM_int/nrow(data_eval)
}

save(log_score, file = "log_score_cml_bs_extreme_100.Rda")


## Boxplot of the out-of-sample censored log scores:
load("log_score_cml_bs_extreme_100.Rda")
log_score_dat <- data.frame(log_score = c(log_score[,1], log_score[,2], 
                                          log_score[,3], log_score[,4],
                                          log_score[,5]),
                             model = c(rep("Cox", 100), rep("Cox_int", 100),
                                       rep("Cforest", 100), rep("CTM", 100),
                                       rep("CTM_int", 100)))
boxplot(log_score ~ model, data = log_score_dat)


