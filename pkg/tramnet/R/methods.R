# logLik method for class "tramnet"

logLik.tramnet <- function(object,
                           parm = coef(object, tol = 0, with_baseline = TRUE),
                           w = NULL, newdata, ...) {
  if (length(list(...)) > 0L)
    warning("additional arguments ignored")

  ctmobj <- .tramnet2ctm(object)
  if (missing(newdata))
    newdata <- .get_tramnet_data(object)
  mltobj <- mlt(ctmobj, data = newdata, dofit = FALSE)
  if (object$model$negative)
    parm <- .flip_sign(object, parm)
  ret <- logLik(mltobj, parm = parm, w = w, ...)
  pen <- .elnet(object)
  attr(ret, "df") <- NA
  class(ret) <- "logLik"
  return(ret - pen)
}

# coef method for class "tramnet"

coef.tramnet <- function(object, with_baseline = FALSE, tol = 1e-6, ...) {
  if (length(list(...)) > 0L)
    warning("additional arguments ignored")
  beta <- c(object$beta[abs(object$beta) > tol])
  if (all(object$x == 0))
    beta <- NULL
  theta <- c(object$theta)
  names(theta) <- names(coef(as.mlt(object$model)))
  if (!with_baseline)
    return(beta)
  return(c(theta, beta))
}

# coef method for class "tramnet_Lm"

coef.tramnet_Lm <- function(object, with_baseline = FALSE, tol = 1e-6,
                            as.lm = FALSE, ...) {
  class(object) <- class(object)[-1L]
  if (!as.lm)
    return(coef(object, with_baseline = with_baseline, tol = tol, ...))
  if (!is.null(object$stratacoef))
    stop("Cannot compute scaled coefficients with strata variables present")
  cf <- coef(object, with_baseline = TRUE, tol = 0, ...)
  cfx <- coef(object, with_baseline = FALSE, tol = 0, ...)
  cfs <- cf[!(names(cf) %in% names(cfx))]
  sd <- 1/cfs[names(cfs) != "(Intercept)"]
  ret <- c(-cf["(Intercept)"], cfx) * sd
  attr(ret, "scale") <- sd
  return(ret)
}


# predict method for class "tramnet"

predict.tramnet <- function(object, newdata = .get_tramnet_data(object), ...) {
  ctmobj <- .tramnet2ctm(object)
  predict(ctmobj, newdata = newdata, ...)
}

# simulate method for class "tramnet"

simulate.tramnet <- function(object, nsim = 1, seed = NULL,
                             newdata = .get_tramnet_data(object),
                             bysim = TRUE, ...) {
  ctmobj <- .tramnet2ctm(object)
  simulate(ctmobj, nsim = nsim, seed = seed,
           newdata = newdata, bysim = bysim, ...)
}

# estfun method for class "tramnet"

estfun.tramnet <- function(object,
                           parm = coef(object, with_baseline = TRUE, tol = 0),
                           w = NULL, newdata, ...) {
  if (any(object$tuning_parm > 0))
    stop("Cannot compute the score for penalised parameters.")
  ctmobj <- .tramnet2ctm(object)
  if (missing(newdata))
    newdata <- .get_tramnet_data(object)
  mltobj <- mlt(ctmobj, data = newdata, dofit = FALSE)
  return(estfun(mltobj, parm = parm, w = w))
}

# residuals method for class "tramnet"

residuals.tramnet <- function(object,
                              parm = coef(object, tol = 0,
                                          with_baseline = TRUE),
                              w = NULL, newdata, ...) {
  ctmobj <- .tramnet2ctm(object)
  if (missing(newdata))
    newdata <- .get_tramnet_data(object)
  mltobj <- mlt(ctmobj, data = newdata, dofit = FALSE)
  return(residuals(object = mltobj, parm = parm,
                   w = w, newdata = newdata, ...))
}

# plot method for class "tramnet"

plot.tramnet <- function(x, newdata,
                         type = c("distribution", "survivor", "density",
                                  "logdensity", "hazard", "loghazard",
                                  "cumhazard", "quantile", "trafo"),
                         q = NULL, prob = 1:(K - 1)/K, K = 50,
                         col = rgb(.1, .1, .1, .1), lty = 1, add = FALSE, ...) {
  ctmobj <- .tramnet2ctm(x)
  if (missing(newdata))
    newdata <- .get_tramnet_data(x)
  plot(ctmobj, newdata = newdata, type = type, q = q, prob = prob, K = K,
       col = col, lty = lty, add = add, ...)
}

# print method for class "tramnet"

print.tramnet <- function(x, ...) {
  print(summary(x, ...))
}

# summary method for class "tramnet"

summary.tramnet <- function(object, ...) {
  tp <- object$model$tram
  if (object$tuning_parm[1] > 0)
    tp <- paste0("Regularized", tp)
  ret <- list(call = getCall(object), convergence = object$result$status,
              type = tp, logLik = logLik(object), coef = coef(object),
              sparsity = .tramnet_sparsity(object),
              tuning_parm = object$tuning_parm)
  class(ret) <- "summary.tramnet"
  ret
}

# print summary method for class "tramnet"

print.summary.tramnet <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nConvergence: ", x$convergence)
  cat("\nType: ", x$type)
  cat("\nLog-Likelihood: ", x$logLik)
  cat("\n")
  cat("\nCoefficients:\n")
  if (!is.null(x$coef)) {
    print(round(x$coef, digits = digits))
  } else {
    print(NULL)
  }
  cat("\nSparsity: ", x$sparsity, "\n")
  cat("\nTuning parameters:\n")
  print(round(x$tuning_parm, digits = digits))
  cat("\n\n")
  invisible(x)
}

# Helper functions

.tramnet2ctm <- function(object) {
  data <- .get_tramnet_data(object)
  cfx <- coef(object, with_baseline = TRUE, tol = 0)
  if (!is.null(object$model$model$model$binteracting)) {
    yBasis <- object$model$model$model$binteracting$iresponse
    iBasis <- object$model$model$model$binteracting$iinteracting
  } else {
    yBasis <- object$model$model$model$bresponse
    iBasis <- NULL
  }

  todistr <- switch(object$model$todistr$name,
                    "minimum extreme value" = "MinExtrVal",
                    "maximum extreme value" = "MaxExtrVal",
                    "normal" = "Normal", "logistic" = "Logistic")
  if (all(object$x == 0)) {
    shifting <- NULL
  } else {
    shifting <-
      as.basis(
        as.formula(
          paste("~", paste(colnames(object$x), collapse = "+"))
        ), data = data, remove_intercept = TRUE
      )
  }
  mod <- ctm(response = yBasis, shifting = shifting,
             interacting = iBasis, todistr = todistr, data = data)
  coef(mod) <- cfx
  return(mod)
}

.get_tramnet_data <- function(object) {
  return(cbind(object$model$data, object$x))
}

.tramnet_sparsity <- function(object, ...) {
  n_coef <- length(coef(object, tol = 0))
  non_zero_coef <- length(coef(object, ...))
  ret <- paste(n_coef, "regression coefficients,", non_zero_coef,
               "of which are non-zero")
  return(ret)
}

.elnet <- function(object) {
  lambda <- object$tuning_parm[1]
  alpha <- object$tuning_parm[2]
  cfx <- coef(object, tol = 0)
  if (is.null(cfx))
    return(0)
  L1 <- sum(abs(cfx))
  L2 <- sqrt(sum(cfx^2))
  ret <- lambda * (0.5 * (1 - alpha) * L2^2 + alpha * L1)
  names(ret) <- NULL
  attr(ret, "L1") <- L1
  attr(ret, "L2") <- L2
  return(ret)
}

.flip_sign <- function(object, parm) {
  stopifnot(inherits(object, "tramnet"))
  wbl <- length(coef(object, tol = 0, with_baseline = TRUE))
  wobl <- length(coef(object, tol = 0))
  idx <- wbl - wobl + seq_len(wobl)
  parm[idx] <- - parm[idx]
  return(parm)
}
