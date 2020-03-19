# Cross validation for "tramnet" models

cvl_tramnet <- function(object, fold = 2, lambda = 0, alpha = 0, folds = NULL,
                        fit_opt = FALSE) {
  if (!inherits(object, "tramnet"))
    stop("Cross validation only for models of class tramnet.")
  df <- .get_tramnet_data(object)
  rsp <- variable.names(object$model, "response")
  n <- nrow(df)
  val_grid <- expand.grid(lambda = lambda, alpha = alpha)
  if (is.null(folds) & !is.null(fold)) {
    folds <- sample(rep(1:fold, ceiling(n/fold)), n)
  } else {
    folds <- round(folds)
    fold <- max(folds)
  }
  out <- .cvl_helper(val_grid = val_grid, df = df, rsp = rsp, folds = folds,
                     object = object, fold = fold, n = n)
  raw_ll <- .get_logLik(out, fold = fold)
  ll_tab <- cbind(val_grid, raw_ll)
  raw_cfx <- .get_cfx(out, fold = fold)
  cfx_tab <- lapply(raw_cfx, function(x) cbind(val_grid, x))
  optimal <- ll_tab[which.max(raw_ll$sum_logLik), ]
  ret <- list(
    logLik_tab = ll_tab,
    optimal = optimal,
    coefficients = cfx_tab,
    folds = folds
  )
  if (fit_opt) {
    optmod <- update(object, lambda = optimal$lambda, alpha = optimal$alpha)
    ret$full_fit <- optmod
  }
  class(ret) <- "cvl_tramnet"
  return(ret)
}

# Helper functions

.cvl_helper <- function(val_grid, df, rsp, folds, object, fold, n) {
  ret <- apply(val_grid, 1, function(pars) {
    message("Performing ", fold, "-fold cross validation")
    sapply(seq_len(fold), function(x, lmb = pars[1], alp = pars[2]) {
      message("Fold: ", x)
      idx <- which(x == folds)
      trn <- df[-idx, , drop = FALSE]
      tst <- df[idx, , drop = FALSE]
      xtrn <- object$x[-idx, , drop = FALSE]
      xtst <- object$x[idx, , drop = FALSE]
      fit <- update.default(object$model, data = trn)
      trnt <- try(tramnet(model = fit, x = xtrn, lambda = lmb, alpha = alp,
                          check_dcp = FALSE, solver = "ECOS"))
      ncfx <- length(coef(object, tol = 0))
      if (inherits(trnt, "try-error")) {
        list(ll = NA,
             cfx = rep(NA, ncfx))
      } else {
        list(ll = logLik(trnt, newdata = tst),
             cfx = coef(trnt, tol = 0))
      }
    }, simplify = FALSE)
  })
  return(ret)
}

.get_logLik <- function(cvl_helper_out, fold) {
  tmp <- sapply(seq_len(fold), function(x) {
    lapply(cvl_helper_out, function(y) y[[x]]$ll)
  }, simplify = FALSE)
  tmp <- lapply(tmp, function(x) do.call(rbind, x))
  ret <- do.call(cbind, tmp)
  colnames(ret) <- paste("logLik_fold", seq_len(fold), sep = "_")
  ret <- data.frame(ret)
  ret$sum_logLik <- apply(ret, 1, sum, na.rm = FALSE)
  return(ret)
}

.get_cfx <- function(cvl_helper_out, fold) {
  tmp <- sapply(seq_len(fold), function(y) {
    lapply(cvl_helper_out, function(x) x[[y]]$cfx)
  }, simplify = FALSE)
  ret <- lapply(tmp, function(x) do.call(rbind, x))
  names(ret) <- paste("coef_fold", seq_len(fold), sep = "_")
  return(ret)
}

.plot_cvl <- function(object, ...) {
  stopifnot(inherits(object, "cvl_tramnet"))
  ll <- object[["logLik_tab"]]
  cfx <- object[["coefficients"]]
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mar = c(5, 5, 2, 2))
  plot(ll$lambda, ll$sum_logLik, ylab = "CV logLik",
       xlab = expression(lambda))
  lapply(cfx, function(cfxx) {
    matplot(x = cfxx[,1], cfxx[,-(1:2)], type = "l",
            xlab = expression(lambda),
            ylab = expression(hat(beta)[j](lambda)), ...)
  })
}
