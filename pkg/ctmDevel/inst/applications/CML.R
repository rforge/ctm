
##********************************************************************
## Chronic myelogenous leukaemia
##********************************************************************

## The following R Code allows to retrace the analysis of the chronic myelogenous
## leukemia data in Section 4 of the main paper. Since we are not able to deliver
## the original data set, we simulated a comparable data set. The relevant
## explanatory variables and the survival times are sampled randomly and
## independently. The number of observations n = 507 equals the number of
## observations in the original data set. The data set is analysed in terms of a
## conditional transformation model and a frailty Cox model. The censored 
## log score is used for model evaluation in a bootstrap approach.

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



#######################################################
### code chunk number 4: Chronic myelogenous leukaemia
#######################################################

library("ctmDevel")
library("prodlim")
library("survival")

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
cml_test <- data.frame(treatment = treatment, riskgroup = riskgroup, 
                       sex = sex, age = age, center = center, 
                       time = sT, status = I(T <= C))



## ************ Conditional Transformation Model ********************


## Define treatment - riskgroup interaction
cml_test$riskgroup_treatment <- interaction(cml_test$treatment, 
                                            cml_test$riskgroup)
ngrid <- sort(unique(cml_test$time))

## Calculate inverse probability of censoring weights
## Calculate Kaplan-Meier estimator of the censoring distribution
G <- prodlim(Surv(time, status) ~ 1, data = cml_test, 
             reverse = TRUE)
ipcw <- inv_weights_graf(ngrid = ngrid, sT = cml_test$time, 
                         event = cml_test$status, G = G)

## Set up model forumla:
## Separate smooth functions over time are estimated for each 
## sex and treatment-riskgroup category. A smooth bivariate 
## surface for time and age is estimated.
## Choose linear base-learners bols for categorical and B-Spline 
## base-learners bbs for metric covariates.
fm <- "bbs(time, df = 2.05, constraint = \"increasing\") ~
                   bbs(age, df = 2.05) +
                   bols(sex, intercept = FALSE, df = 2.05) +
                   bols(riskgroup_treatment, intercept = FALSE, 
                        df = 3)"
fm <- as.formula(fm)

## Estimate CTM:
## Include random intercepts (constant over time) for each study 
## center. 
options(mboost_useMatrix = FALSE) 
ctm_cml_int <- mctm(fm, data = cml_test, family = Binomial_gumbel(),
                    constant = "brandom(center, df = 2.05^2)", 
                    ngrid = ngrid, weights = ipcw, control = 
                    boost_control(nu = 0.5, mstop = 2000, 
                    trace = TRUE))


## ***************** Extract estimated survival probabilities
##           (E.g. to plot category-specific survival curves)

## Calculate estimated survival probabilites for each treatment-
## riskgoup combination, each sex and age = 36, 48, 58. Set all 
## other covariates to their reference values: sex = female, 
## treatment-riskgroup = IFN-alpha.0 (IFN-alpha treatment in the 
## low risk group), center = 1 and age = 45.
vn <- names(variable.names(ctm_cml_int))

## Treatment-riskgroup 
length_grid <- length(unique(cml_test$time)) + 1
pred_IFN0 <- pred_BUS0 <- pred_HU0 <- 
matrix(0, nrow = length_grid, ncol = length(vn))
pred_IFN1 <- pred_BUS1 <- pred_HU1 <- 
matrix(0, nrow = length_grid, ncol = length(vn))
pred_IFN2 <- pred_BUS2 <- pred_HU2 <- 
matrix(0, nrow = length_grid, ncol = length(vn))
for(i in 1:length(vn)){
  pred_IFN0[,i] <- predict(ctm_cml_int, newdata = data.frame(
                          riskgroup_treatment = factor("IFN-alpha.0", 
                          levels = levels(cml_test$riskgroup_treatment)),
                          sex = factor(0, levels = c(0,1)), age = 45,
                          center = factor(1, levels = 
                          levels(cml_test$center))), which = vn[i], 
                          anno = TRUE,
                          y = c(0, sort(unique(cml_test$time))))$p
  pred_BUS0[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("BUS.0", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(0, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
  pred_HU0[,i] <- predict(ctm_cml_int, newdata = data.frame(
                          riskgroup_treatment = factor("HU.0", 
                          levels = levels(cml_test$riskgroup_treatment)),
                          sex = factor(0, levels = c(0,1)), age = 45,
                          center = factor(1, levels = 
                          levels(cml_test$center))),
                          which = vn[i], anno = TRUE,
                          y = c(0, sort(unique(cml_test$time))))$p
  pred_IFN1[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("IFN-alpha.1", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(0, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
  pred_BUS1[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("BUS.1", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(0, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
  pred_HU1[,i] <- predict(ctm_cml_int, newdata = data.frame(
                          riskgroup_treatment = factor("HU.1", 
                          levels = levels(cml_test$riskgroup_treatment)),
                          sex = factor(0, levels = c(0,1)), age = 45,
                          center = factor(1, levels = 
                          levels(cml_test$center))),
                          which = vn[i], anno = TRUE,
                          y = c(0, sort(unique(cml_test$time))))$p
  pred_IFN2[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("IFN-alpha.2", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(0, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
  pred_BUS2[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("BUS.2", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(0, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
  pred_HU2[,i] <- predict(ctm_cml_int, newdata = data.frame(
                          riskgroup_treatment = factor("HU.2", 
                          levels = levels(cml_test$riskgroup_treatment)),
                          sex = factor(0, levels = c(0,1)), age = 45,
                          center = factor(1, levels = 
                          levels(cml_test$center))),
                          which = vn[i], anno = TRUE,
                          y = c(0, sort(unique(cml_test$time))))$p
}
## Estimated survival probabilities for treatment-riskgroup categories
## treatment = IFN-alpha, risk = low
pis_IFN0 <- pgumbel(rowSums(pred_IFN0)) 
## treatment = BUS, risk = low
pis_BUS0 <- pgumbel(rowSums(pred_BUS0)) 
## treatment = HU, risk = low
pis_HU0 <- pgumbel(rowSums(pred_HU0))
## treatment = IFN-alpha, risk = mediate   
pis_IFN1 <- pgumbel(rowSums(pred_IFN1)) 
## treatment = BUS, risk = mediate
pis_BUS1 <- pgumbel(rowSums(pred_BUS1))
## treatment = HU, risk = mediate 
pis_HU1 <- pgumbel(rowSums(pred_HU1))  
## treatment = IFN-alpha, risk = high 
pis_IFN2 <- pgumbel(rowSums(pred_IFN2))
## treatment = BUS, risk = high 
pis_BUS2 <- pgumbel(rowSums(pred_BUS2))
## treatment = HU, risk = high 
pis_HU2 <- pgumbel(rowSums(pred_HU2))   

## Sex
pred_sex0 <- pred_sex1 <- 
matrix(0, nrow = length(sort(unique(cml_test$time)))+1, ncol = length(vn))
for(i in 1:length(vn)){
  pred_sex0[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("IFN-alpha.0", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(0, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
  pred_sex1[,i] <- predict(ctm_cml_int, newdata = data.frame(
                           riskgroup_treatment = factor("IFN-alpha.0", 
                           levels = levels(cml_test$riskgroup_treatment)),
                           sex = factor(1, levels = c(0,1)), age = 45,
                           center = factor(1, levels = 
                           levels(cml_test$center))),
                           which = vn[i], anno = TRUE,
                           y = c(0, sort(unique(cml_test$time))))$p
}
## Estimated survival probabilites
pis_sex0 <- pgumbel(rowSums(pred_sex0)) ## female
pis_sex1 <- pgumbel(rowSums(pred_sex1)) ## male

## Age
pred_age36 <- pred_age48 <- pred_age58 <-
matrix(0, nrow = length(sort(unique(cml_test$time))) +1, ncol = length(vn))
for(i in 1:length(vn)){
  pred_age36[,i] <- predict(ctm_cml_int, newdata = data.frame(
                            riskgroup_treatment = factor("IFN-alpha.0", 
                            levels = levels(cml_test$riskgroup_treatment)),
                            sex = factor(0, levels = c(0,1)), age = 36,
                            center = factor(1, levels = 
                            levels(cml_test$center))),
                            which = vn[i], anno = TRUE,
                            y = c(0, sort(unique(cml_test$time))))$p
  pred_age48[,i] <- predict(ctm_cml_int, newdata = data.frame(
                            riskgroup_treatment = factor("IFN-alpha.0", 
                            levels = levels(cml_test$riskgroup_treatment)),
                            sex = factor(0, levels = c(0,1)), age = 48,
                            center = factor(1, levels = 
                            levels(cml_test$center))),
                            which = vn[i], anno = TRUE,
                            y = c(0, sort(unique(cml_test$time))))$p
  pred_age58[,i] <- predict(ctm_cml_int, newdata = data.frame(
                            riskgroup_treatment = factor("IFN-alpha.0", 
                            levels = levels(cml_test$riskgroup_treatment)),
                            sex = factor(0, levels = c(0,1)), age = 58,
                            center = factor(1, levels = 
                            levels(cml_test$center))),
                            which = vn[i], anno = TRUE,
                            y = c(0, sort(unique(cml_test$time))))$p
}
## Estimated survival probabilites
pis_age36 <- pgumbel(rowSums(pred_age36))
pis_age48 <- pgumbel(rowSums(pred_age48))
pis_age58 <- pgumbel(rowSums(pred_age58))



## ******************************* Frailty Cox model ******************* 


cml_cox_flex <- coxph(Surv(time, status) ~ treatment + sex + age + 
                           riskgroup + treatment:riskgroup + 
                           frailty.gaussian(center, sparse = FALSE),
                           data = cml_test)
summary(cml_cox_flex)

##******************* Corresponding estimated survivor functions
## The extraction of the survival probabilites needs to be done by 
## hand because the survfit-function does not work properly since no 
## estimate for center = 1 is delivered.

## Treatment-Riskgroup 
## fixed age = 45, sex = 0, center = 1
newdat <- data.frame(age = 45, sex = 0, center = factor(1, levels = 
                     levels(cml_test$center)), riskgroup = 
                     factor(c(0,0,0,1,1,1,2,2,2), levels = c(0,1,2)),
                     treatment = factor(c("IFN_alpha", "BUS", "HU", 
                     "IFN_alpha","BUS","HU","IFN_alpha", "BUS", "HU"),
                     levels = c("IFN_alpha","BUS", "HU")))
## Set up model matrix
modmat <- model.matrix(~treatment+sex+age+riskgroup+treatment:riskgroup
                        + center, newdat)[,-1] 
## add reference category center1 to the model matrix
center1 <- rep(0, nrow(modmat))
center1[which(newdat$center == 1)] <- 1
p <- 10 + length(levels(cml_test$center)) - 1
modelmat <- cbind(modmat[,1:6], center1, modmat[,7:p])
## Calculate the linear predictor
lp <- modelmat %*% coef(cml_cox_flex)
## fitted hazard ratios
HR <- exp(lp)
## baseline hazard
H <- basehaz(cml_cox_flex, centered = FALSE)
hazards <- H$hazard
## Add grid point 0
hazards <- c(0, hazards)
timepoints <- c(0, H$time)
## Calculate survival probabilities for each individual at each grid point
surv_probs_tr_risk <- matrix(0, nrow = length(HR), ncol = length(hazards))
for(i in 1:length(HR)){
  surv_probs_tr_risk[i,] <- exp(-hazards*HR[i])}

## Sex
## fixed age = 45, risk = 0, treatment = IFN-alpha, center = 1
newdat <- data.frame(age = 45, sex = factor(c(0,1), levels = c(0,1)),
                     center = factor(1, levels = levels(cml_test$center)),
                     riskgroup = factor(0, levels = c(0,1,2)),
                     treatment = factor("IFN_alpha", levels = c("IFN_alpha","BUS", "HU")))
modmat <- model.matrix(~treatment+sex+age+riskgroup+treatment:riskgroup
                        + center, newdat)[,-1]
center1 <- rep(0, nrow(modmat))
center1[which(newdat$center == 1)] <- 1
p <- 10 + length(levels(cml_test$center)) - 1
modelmat <- cbind(modmat[,1:6], center1, modmat[,7:p])
lp <- modelmat %*% coef(cml_cox_flex)
HR <- exp(lp)
H <- basehaz(cml_cox_flex, centered = FALSE)
hazards <- H$hazard
hazards <- c(0, hazards)
timepoints <- c(0, H$time)
surv_probs_sex <- matrix(0, nrow = length(HR), ncol = length(hazards))
for(i in 1:length(HR)){
    surv_probs_sex[i,] <- exp(-hazards*HR[i])}

## Age
newdat <- data.frame(age = c(36, 48, 58), sex = factor(0, levels = c(0,1)),
                     center = factor(1, levels = levels(cml_test$center)),
                     riskgroup = factor(0, levels = c(0,1,2)),
                     treatment = factor("IFN_alpha", levels = c("IFN_alpha","BUS", "HU")))
modmat <- model.matrix(~treatment+sex+age+riskgroup+treatment:riskgroup
                        + center, newdat)[,-1]
center1 <- rep(0, nrow(modmat)) 
center1[which(newdat$center == 1)] <- 1
p <- 10 + length(levels(cml_test$center)) - 1
modelmat <- cbind(modmat[,1:6], center1, modmat[,7:p])
lp <- modelmat %*% coef(cml_cox_flex)
HR <- exp(lp)
H <- basehaz(cml_cox_flex, centered = FALSE)
hazards <- H$hazard
hazards <- c(0, hazards)
timepoints <- c(0, H$time)
surv_probs_age <- matrix(0, nrow = length(HR), ncol = length(hazards))
for(i in 1:length(HR)){
    surv_probs_age[i,] <- exp(-hazards*HR[i])}




## ************************* Model evaluation ***************************


## Evaluation of the CTM and the frailty Cox model performance for new 
## observations.
## The estimation and the evaluation data set are determined in terms 
## of a bootstrap approach. B = 50 bootstrap replications are conducted 
## and model performance is measured in terms of the censored log-score.

## Generate estimation and evaluation data sets
nboot <- 50  ## number of bootstrap replications
set.seed(24)
indx <- vector(mode = "list", length = nboot)
for(i in 1:nboot){
  ## Step 1: Sample n1 observations with replacement.
  n1 <- nrow(cml_test) - length(levels(cml_test$center))
  indx1 <- sample(1:nrow(cml_test), n1, replace = TRUE)
  ## Step2: Sample one observation from EACH study centre.
  ## This procedure guarantees that all study centres are represented 
  ## in each booststrap data set. 
  level <- unique(cml_test$center)
  indx2 <- rep(0, length(level))
     for(j in 1:length(level)){
     ## Choose all individuals from the specific study centre:
     indx_subs <- which(cml_test$center == level[j])
     if(length(indx_subs) == 1) indx2[j] <- indx_subs
     else{ indx2[j] <- sample(indx_subs,1)}}
  indx[[i]] <- c(indx1, indx2)
}

## Calculate the censored log-score for each out-of-bootstrap sample:
log_score <- matrix(0, nrow = nboot, ncol = 2)
colnames(log_score) <- c("log_score CTM", "log-score Cox")
rownames(log_score) <- 1:nboot
for(i in 1:nboot){
   ## Define estimation and evaluation data set:
   data_est <- cml_test[indx[[i]],]
   data_eval <- cml_test[-unique(indx[[i]]),]

   ##************************* Estimate CTM
   ## Calculate the evaluation grid points:  
   ngrid_eval <- sort(unique(data_eval$time))
   ## Calculate IPCWs:
   G <- prodlim(Surv(time, status) ~ 1, data = data_est, reverse = TRUE)
   ipcw_eval <- inv_weights_graf(ngrid = ngrid_eval, sT = data_est$time,
                                 event = data_est$status, G = G)
   ipcw_mat <- matrix(ipcw_eval, ncol = length(ngrid_eval))

   ## Estimate CTM for the estimation data set:
   data_est$riskgroup_treatment <- interaction(data_est$treatment, 
                                               data_est$riskgroup)
   fm <- "bbs(time, df = 2.05, constraint = \"increasing\") ~ 
                   bbs(age, df = 2.05) + 
                   bols(sex, intercept = FALSE, df = 2.05) +
                   bols(riskgroup_treatment, intercept = FALSE, df = 3)"
   fm <- as.formula(fm)
   ctm_cml <- mctm(fm, data = data_est, family = Binomial_gumbel(),
                   constant = "brandom(center, df = 2.05^2)",
                   ngrid = ngrid_eval, weights = ipcw_eval,
                   control = boost_control(nu = 0.5, mstop = 2000,
                   trace = FALSE))
   ## Estimated survival probabilites:
   pis_eval <- predict(ctm_cml, newdata = data_eval,
                       y = sort(unique(data_eval$time)), 
                       type = "response")
   ## NOTE: The censored log-score is only evaluated at the event and 
   ## censoring time points of the evaluation data set!
   pis_eval_mat <- matrix(pis_eval, ncol = 
                          length(sort(unique(data_eval$time))))
   ## Calculate censored log_score:
   log_score_temp <- rep(0, nrow(data_eval))
   ipcw_eval_mat <- matrix(ipcw_eval, ncol = length(ngrid_eval))

   for(l in 1:nrow(data_eval)){
       log_score_temp[l] <- log_score_cens(data_eval$time[l],
                                         sort(unique(data_eval$time)),
                                         pis_eval_mat[l,], ipcw_eval_mat[l,])
   }
   log_score[i,1] <- sum(log_score_temp)/nrow(data_eval)

   ## ***************************** Estimate frailty Cox model
   cml_cox_flex <- coxph(Surv(time, status) ~ treatment + sex + age + 
                              riskgroup + treatment:riskgroup + 
                              frailty.gaussian(center, sparse = FALSE), 
                              data = data_est)
   ## Extract corresponding survival probabilities:
   ## Calculate linear predictor:
   modmat <- model.matrix(~treatment+sex+age+riskgroup+treatment:riskgroup
                           + center, data_eval)[,-1]
   center1 <- rep(0, nrow(modmat))
   indx_center1 <- which(data_eval$center == 1)
   center1[indx_center1] <- 1
   p <- 10 + length(levels(data_est$center)) -1
   modelmat <- cbind(modmat[,1:6], center1, modmat[,7:p])
   lp <- modelmat %*% coef(cml_cox_flex)
   ## fitted hazard ratios
   HR <- exp(lp)
   ## baseline hazard
   H <- basehaz(cml_cox_flex, centered = FALSE)
   h0 <- rep(0, length = length(sort(unique(data_eval$time))))
   for(idx in 1:length(h0)){
       if(length(which(H$time <= ngrid_eval[idx])) == 0) k <- 1
       else{ k <- max(which(H$time <= ngrid_eval[idx]))}
       h0[idx] <- H$hazard[k]
   }
   pis_cox <- matrix(0, nrow = length(HR), ncol = length(h0))
   for(idx in 1:length(HR)){
       pis_cox[idx,] <- exp(-h0*HR[idx])
   }
   pis_cox <- 1 - pis_cox ## We need probabilities of the distribution 
                          ## function.
   ## Calculate censored log score:
   log_score_temp_cox <- c()
   for(l in 1:nrow(data_eval)){
       log_score_temp_cox[l] <- log_score_cens(data_eval$time[l],
                                            sort(unique(data_eval$time)),
                                            pis_cox[l,], ipcw_eval_mat[l,])
   }
   log_score[i,2] <- sum(log_score_temp_cox)/nrow(data_eval)
}
## Box-plot of the censored log scores:
dat <- data.frame(log_score = c(log_score[,1], log_score[,2]),
                  model = c(rep("CTM", nrow(log_score)), rep("Cox", nrow(log_score))))
boxplot(log_score ~ model, data = dat)


