
library("tbm")
library("tram")
library("trtf")
library("TH.data")

set.seed(290875)

NCORES <- 10
B <- 100
M <- 1000

FUN <- function(model, fm, ldata, control, method = quote(mboost::mboost), 
                folds) {

  w <- weights(model)
  fmin <- fm
  ### this is where model.frame.default in blackboost looks for weights
  if (length(fm) ==  2) {
       fm <- fmin[["ctm"]]
  } else {
       assign("w", w, environment(fm))
  }
  ctm <- ctmboost(model, formula = fm, data = ldata, weights = w, 
                      control = bctrl, method = method)
  cv_ctm <- cvrisk(ctm, folds = folds, mc.cores = NCORES)
  ctm <- ctm[mstop(cv_ctm)]

  if (length(fmin) == 2) fm <- fmin[["tram"]]
  tram <- stmboost(model, formula = fm, data = ldata, weights = w, 
                    control = control, method = method)
  cv_tram <- cvrisk(tram, folds = folds, mc.cores = NCORES)
  tram <- tram[mstop(cv_tram)]

  risk <- cbind(ctm = -cv_ctm[,mstop(cv_ctm) + 1],
                tram = -cv_tram[,mstop(cv_tram) + 1])

  list(ctm = ctm, tram = tram, risk = risk, 
       cv_ctm = cv_ctm, cv_tram = cv_tram)
}

FUN2 <- function(model, fm, ldata, tcontrol, fcontrol, folds) {

  ret <- matrix(NA, ncol = 2, nrow = ncol(folds))
  colnames(ret) <- c("tree", "forest")

  for (j in 1:ncol(folds)) {
      print(j)
      ll <- ldata[rep(1:nrow(ldata), fd[,j]),]
      tt <- ldata[fd[,j] == 0,]
      tr_theta <- trafotree(model, formula = fm, data = ll, 
                            control = tcontrol)
      ret[j, "tree"] <- logLik(tr_theta, newdata = tt) / nrow(tt)

      tf_theta <- try(traforest(model, formula = fm, data = ll, 
                            control = fcontrol, ntree = 100, trace = TRUE,
                            cores = NCORES))
      if (inherits(tf_theta, "try-error")) next()
      cf <- try(predict(tf_theta, newdata = tt, type = "coef", cores = NCORES))
      if (inherits(cf, "try-error")) next()
      tmp <- do.call("rbind", cf)
      idx <- which(abs(tmp) > 10, arr.ind = TRUE)[,1]
      if (length(idx) > 1) 
          cf[idx] <- predict(tf_theta, newdata = tt[idx,,drop =FALSE], type = "coef")
      ret[j, "forest"] <- logLik(tf_theta, newdata = tt, coef = cf) / nrow(tt)
  }

  ret
}
