library("cotram") 
library("MASS")
library("multcomp")


## set-up plots
library("lattice")
library("gridExtra")
library("colorspace")
col <- qualitative_hcl(3, c(240, 0), l = 60)


## loading data
if (!file.exists("analysis/DVC.rda")) {
  download.file("https://zenodo.org/record/17179/files/DVC.tgz", "DVC.tgz")
  untar("DVC.tgz", file = "analysis/DVC.rda")
}
load("analysis/DVC.rda")


## settings
set.seed(2925)
loc <- Sys.setlocale("LC_ALL", "en_US.UTF-8")
rm(loc)


## data frame
df <- data.frame(wild = as.vector(obs[,,"wild"]),
                 day = seq(from = start, to = end,by = "30 min")[1:prod(dim(obs)[1:2])])
df$weekday <- factor(format(df$day, "%A"),
                     levels = c("Monday", "Tuesday", "Wednesday",
                                "Thursday", "Friday", "Saturday", "Sunday"))
df$time <- as.numeric(difftime(df$day, start, unit = "days"))
df$year <- factor(format(df$day, "%Y"))
df$Datum <- factor(format(df$day, "%Y-%m-%d"))
df$daytime <- cut(as.numeric(difftime(df$day, trunc(df$day, "day"), unit = "hours")), 
                  breaks = c(-Inf, 12, Inf), labels = c("am", "pm"))


## extract sunrise, sunset
w <- weekdays[, c("Datum", "SAofficial", "SUofficial", "ArbZeitFaktor")]
w$Datum <- factor(format(w$Datum, "%Y-%m-%d"))

df <- merge(df, w, by = "Datum")
df$weekday[df$ArbZeitFaktor == 0] <- "Sunday"

a <- as.numeric(with(df, difftime(day, SAofficial, unit = "mins")))
df$sunrise <- cut(a, breaks = c(-Inf, -120, -15, 120, Inf), 
                  labels = c("night", "pre.sunrise", "post.sunrise", "day")) 

a <- as.numeric(with(df, difftime(day, SUofficial, unit = "mins")))
df$sunset <- cut(a, breaks = c(-Inf, -120, -15, 120, Inf), 
                 labels = c("day", "pre.sunset", "post.sunset", "night")) 

df$hours <- with(df, interaction(sunrise, sunset))[, drop = TRUE]
levels(df$hours) <- c("night", "pre.sunrise", "post.sunrise", "day",
                      "pre.sunset", "post.sunset", "night")
df$hours <- relevel(df$hours, "day")
df$daytime <- with(df, interaction(hours, daytime)[, drop = TRUE])


## harmonics
sfm <- function(timevar, freq = 365, S = 10) {
  S <- 1:S * 2
  paste("I(", rep(c("sin", "cos"), length(S)), "(",
        rep(S, rep(2, length(S))), "* pi *", timevar, "/", freq, "))",
        collapse = "+")
}

Xtime <- model.matrix(as.formula(paste("~", sfm("time"))), data = df)[,-1]
Xtime <- as.data.frame(Xtime)
colnames(Xtime) <- paste0("tvar", 1:ncol(Xtime))
df <- cbind(df, Xtime)


## training and validation data-set
trainID <-  as.character(2002:2009)
df_train <- subset(df, year %in% trainID)
df_train$year <- droplevels(df_train$year)
df_test <- subset(df, !year %in% trainID)


## last value carried forward
df_test$year <- "2009"
df_test$year <- factor(df_test$year, levels = levels(df_train$year))


## model formula
rhs <- paste("year + daytime * weekday +", 
             paste("daytime:", colnames(Xtime), "", collapse = "+"))

fm <- as.formula(paste("wild ~", rhs))


## Poisson GLM
mp <- glm(fm, data = df_train, family = poisson(link = "log"))

### Poisson GLM: out-of-sample log-likelihood
sum(dpois(df_test$wild,
          lambda = predict(mp, newdata = df_test, type = "response"),
          log = TRUE))


## neg. bin model
mnb <- glm.nb(fm, data = df_train, link = "log")

### neg. bin model: out-of-sample log-likelihood
sum(dnbinom(df_test$wild,
            mu = predict(mnb, newdata = df_test, type = "response"),
            size = mnb$theta, log = TRUE))


### Cox count transformation model
mcll <- cotram(fm, data = df_train, method = "cloglog", prob = .99)

### Cox count transformation model: out-of-sample log-likelihood
logLik(mcll, newdata = df_test)


## estimated risk and 95% simultaneous confidence bands
ci.FUN <- function(model){
  ret <- vector(mode = "list", length = nlevels(df$daytime))
  names(ret) <- levels(df$daytime)
  
  for (d in levels(df$daytime)) {
    nd_am <- subset(df_train, day > as.POSIXlt("2009-01-01 00:00:00", tz = "UTC") &
                      day < as.POSIXlt("2009-12-31 24:00:00", tz = "UTC")) 
    nd_am <- subset(nd_am, daytime == d)
    nd_am <- subset(nd_am, !duplicated(Datum))
    nd_am$weekday <- factor("Monday",
                            levels = c("Monday", "Tuesday", "Wednesday",
                                       "Thursday", "Friday", "Saturday","Sunday"))
    
    ci_d.FUN <- function(model, nd, FUN = exp, vcov) {
      
      K <- model.matrix(fm, data = nd[yrw <- seq(from = 1, to = 365, by = 15),])
      K[, c("year2009", colnames(K)[grep("(Intercept)", colnames(K))])] <- 0
      if("cotram" %in% class(model)) K <- K[,-1] # remove intercept for cotram
      ci <- confint(glht(model, linfct = K, vcov. = vcov))
      alpha <- attr(ci$confint, "calpha")
      
      K <- model.matrix(fm, data = nd)
      if("cotram" %in% class(model)) K <- K[,-1] # remove intercept for cotram
      ci <- confint(glht(model, linfct = K, vcov. = vcov), calpha = alpha)
      
      eJan01 <- ci$confint[1, "Estimate"]
      ci$confint[,1] <- ci$confint[,1] - eJan01
      ci$confint[,2] <- ci$confint[,2] - eJan01
      ci$confint[,3] <- ci$confint[,3] - eJan01
      data.frame(FUN(ci$confint), day = nd$day)
    }
    
    ret[[d]] <- ci_d.FUN(model = model, nd = nd_am, vcov = vcov(model))
    ret[[d]]$model <- NULL
  }
  for (d in names(ret)) ret[[d]]$daytime <- d
  return(ret)
}

fit_mp <- ci.FUN(model = mp)
fit_mnb <- ci.FUN(model = mnb)
fit_mcll <- ci.FUN(model = mcll)

## plots of estimated risk and 95% simultaneous confidence bands
plotFUN <- function(fit){
  tmp <- do.call("rbind", fit)
  tmp$daytime <- ordered(as.character(tmp$daytime),
                         levels = c("night.am", "pre.sunrise.am", "post.sunrise.am",
                                    "day.am", "day.pm", "pre.sunset.pm",
                                    "post.sunset.pm", "night.pm"),
                         labels = c("Night (am)", "Pre-sunrise", "Post-sunrise",
                                    "Day (am)", "Day (pm)", "Pre-sunset",
                                    "Post-sunset", "Night (pm)"))
  return(tmp)
}

pp <- xyplot(cbind(lwr, Estimate, upr) ~ day | daytime, data = plotFUN(fit = fit_mp),
             ylab = "Expectation ratio", xlab = "Day of year", main = "(A) Poisson model",
             panel = function(x, y, ...) {
               panel.abline(h = 1, col = "gray70")
               col <- rgb(.3, .3, .3, .3)
               y <- matrix(y, ncol = 3)
               panel.xyplot(x, y[,2], type = "l", col = "black")
               panel.polygon(c(x, rev(x)), c(y[, 1], rev(y[, 3])),
                             col = col, border = NA)},
             layout = c(4, 2), ylim = c(0, 5.5))


pnb <- xyplot(cbind(lwr, Estimate, upr) ~ day | daytime, data = plotFUN(fit = fit_mnb),
              ylab = "Expectation ratio", xlab = "Day of year", main = "(B) Negative binomial model",
              panel = function(x, y, ...) {
                panel.abline(h = 1, col = "gray70")
                col <- rgb(.3, .3, .3, .3)
                y <- matrix(y, ncol = 3)
                panel.xyplot(x, y[,2], type = "l", col = "black")
                panel.polygon(c(x, rev(x)), c(y[, 1], rev(y[, 3])),
                              col = col, border = NA)},
              layout = c(4, 2), ylim = c(0, 5.5))


pcll <- xyplot(cbind(lwr, Estimate, upr) ~ day | daytime, data = plotFUN(fit = fit_mcll),
               ylab = "Hazard ratio", xlab = "Day of year", main = "(C) Cox count transformation model",
               panel = function(x, y, ...) {
                 panel.abline(h = 1, col = "gray70")
                 col <- rgb(.3, .3, .3, .3)
                 y <- matrix(y, ncol = 3)
                 panel.xyplot(x, y[,2], type = "l", col = "black")
                 panel.polygon(c(x, rev(x)), c(y[, 1], rev(y[, 3])),
                               col = col, border = NA)},
               layout = c(4, 2), ylim = c(0, 5.5))

grid.arrange(pp, pnb, pcll, nrow = 3)


## conditional distribution functions
cdf.FUN <- function(d, q = 0:max(df_train$wild)){
  
  nd <- subset(df, day == as.POSIXlt(d, tz = "UTC"))
  nd$year <- factor(nd$year, levels = levels(df_train$year))
  nd$wild <- NULL
  
  ret <- data.frame(q = q)
  
  ret$p <- ppois(q, lambda = predict(mp, newdata = nd, type = "response"))
  ret$nb <- pnbinom(q, mu = predict(mnb, newdata = nd, type = "response"),
                    size = mnb$theta)
  ret$cll <- c(predict(mcll, newdata = nd, q = q, type = "distribution"))
  
  ret$start <- nd$day
  ret$end <- df[which(df$day == as.POSIXlt(d, tz = "UTC")) + 1, "day"]
  ret$wild <- df[which(df$day == as.POSIXlt(d, tz = "UTC")), "wild"]
  
  return(ret)
}

cdf <- lapply(c("2009-08-29 20:00:00", "2009-04-21 04:30:00",
                "2009-01-01 01:00:00", "2009-02-25 18:00:00"),
              cdf.FUN)

## plots of conditional distribution functions
tmp <- do.call("rbind", cdf)

tmp$interval <- paste(tmp$start, "-", format(tmp$end, "%H:%M:%S"), "(UTC+1)")
tmp$interval <- ordered(as.character(tmp$interval), levels = unique(tmp$interval))

xyplot(cbind(p, nb, cll) ~ q | interval, data = tmp,
       xlab = "Deer-vehicle collisions", ylab = "Conditional distribution", 
       panel = function(x, y, subscripts, ...) {
         panel.abline(h = c(0, 1), col = "gray70", lty = 3)
         y <- matrix(y, ncol = 3)
         for (i in 1:3){
           panel.xyplot(x, y[,i], type = "s", col = col[i])
           panel.xyplot(x, y[,i], type = "p", pch = 19, col = col[i])
         }
         panel.abline(v = tmp$wild[subscripts], col = "black", lwd = 1, lty = 1)},
       key =  list(columns = 1, x = 0.98, y = 0.54, corner = c(1, 0),
                   lines = list(col = col, type = "l", lwd = 1.3),
                   text = list(lab = c("Poisson",
                                       "Negative binomial", 
                                       "Cox count transformation"))),
       layout = c(2, 2)
)



## simultaneous confidence bands for cdf
d2 <- "2009-04-21 04:30:00"
nd <- model.frame(mcll)[which(df$day == as.POSIXlt(d2, tz = "UTC")),]
cb <- as.data.frame(confband(mcll, newdata = nd, type = "distribution"))

cb$start <- d2
cb$end <- df[which(df$day == as.POSIXlt(d2, tz = "UTC")) + 1, "day"]
cb$wild <- df[which(df$day == as.POSIXlt(d2, tz = "UTC")), "wild"]
cb$interval <- paste(cb$start, "-", format(cb$end, "%H:%M:%S"), "(UTC+1)")


## plots of simultaneous confbands for cdf
xyplot(cbind(lwr, Estimate, upr) ~ q | interval, data = cb,
       xlab = "Deer-vehicle collisions", ylab = "Conditional distribution",
       panel = function(x, y, ...) {
         panel.abline(h = c(0, 1), col = "gray70", lty = 3)
         
         y <- matrix(y, ncol = 3)
         est <- y[, 2]
         lwr <- y[, 1]
         upr <- y[, 3]
         polyx <- c(x[1], rep(x[2:length(x)], each = 2),
                    x[length(x)], rep(rev(x)[2:length(rev(x))], each = 2))
         polyy <-  c(rep(lwr[1:(length(lwr)-1)], each = 2), rev(upr)[1],
                     rep(rev(upr)[2:length(rev(upr))], each = 2), lwr[1])
         
         panel.polygon(polyx, polyy, col = rgb(.3, .3, .3, .3), border = NA)
         panel.xyplot(x, est, type = "s", col = col[3])
         panel.xyplot(x, est, type = "p", pch = 19, col = col[3])
         
         panel.abline(v = cb$wild, col = "black", lwd = 1, lty = 1)},
       ylim = c(0 - .1, 1 + .1))

warnings()
