# Profiling tuning parameters

prof_lambda <- function(model, min_lambda = 0, max_lambda = 15, nprof = 5,
                        as.lm = FALSE) {
  stopifnot(inherits(model, "tramnet"))
  stopifnot(max_lambda > min_lambda)
  stopifnot(nprof > 0)
  if (min_lambda == 0)
    lambdas <- c(0, 10^seq(-1, log10(max_lambda), length.out = nprof)[-1])
  else
    lambdas <- 10^seq(log10(min_lambda), log10(max_lambda), length.out = nprof)

  cfx <- list()
  lls <- list()

  for (idx in seq_along(lambdas)) {
    message("Step ", idx, "/", length(lambdas),
            " at lambda = ", round(lambdas[idx], 2))
    mod <- try(update(model, lambda = lambdas[idx]))
    if (inherits(mod, "try-error")) {
      cfx[[idx]] <- rep(NA, length(coef(model, tol = 0, as.lm = as.lm)))
      lls[[idx]] <- NA
    } else {
      cfs <- coef(mod, with_baseline = FALSE, tol = 0, as.lm = as.lm)
      cfx[[idx]] <- cfs[names(cfs) != "(Intercept)"]
      lls[[idx]] <- logLik(mod)
    }
  }

  ret <- list(
    lambdas = lambdas,
    cfx = do.call(rbind, cfx),
    lls = do.call(c, lls)
  )

  class(ret) = "prof_lambda"

  return(ret)
}

# Profiling tuning parameters

prof_alpha <- function(model, min_alpha = 0, max_alpha = 1, nprof = 5,
                       as.lm = FALSE) {
  stopifnot(inherits(model, "tramnet"))
  stopifnot(.check_bounds(min_alpha, max_alpha))
  stopifnot(nprof > 0)
  alphas <- seq(min_alpha, max_alpha, length.out = nprof)
  cfx <- list()
  lls <- list()

  for (idx in seq_along(alphas)) {
    message("Step ", idx, "/", length(alphas),
            " at alpha = ", round(alphas[idx], 2))
    mod <- try(update(model, alpha = alphas[idx]))
    if (inherits(mod, "try-error")) {
      cfx[[idx]] <- rep(NA, length(coef(model, tol = 0, as.lm = as.lm)))
      lls[[idx]] <- NA
    } else {
      cfs <- coef(mod, with_baseline = FALSE, tol = 0, as.lm = as.lm)
      cfx[[idx]] <- cfs[names(cfs) != "(Intercept)"]
      lls[[idx]] <- logLik(mod)
    }
  }

  ret <- list(
    alphas = alphas,
    cfx = do.call(rbind, cfx),
    lls = do.call(c, lls)
  )

  class(ret) = "prof_alpha"

  return(ret)
}

# Plot profiles for "profile_*_tramnet" classes

plot_path <- function(object, plot_logLik = FALSE, ...) {
  if (inherits(object, "prof_lambda"))
    .plot_lpath(object, plot_logLik, ...)
  else
    if (inherits(object, "prof_alpha"))
      .plot_apath(object, plot_logLik, ...)
  else
    stop("plot_path() needs an object of class prof_lambda or prof_alpha")
}

.plot_lpath <- function(object, plot_logLik, ...) {
  if (plot_logLik)
    plot(object$lambdas, object$lls, xlab = expression(lambda),
         ylab = "logLik", ...)
  opar <- par("mar")
  on.exit(par(opar))
  par(mar = c(5, 5, 2, 2))
  matplot(object$lambdas, object$cfx, type = "l", xlab = expression(lambda),
          ylab = expression(hat(beta)[j](lambda)), ...)
  text(x = min(object$lambdas) + 0.05 * abs(diff(range(object$lambdas))),
       y = object$cfx[which.min(object$lambdas), ],
       labels = colnames(object$cfx), ...)
}

.plot_apath <- function(object, plot_logLik, ...) {
  if (plot_logLik)
    plot(object$alphas, object$lls, xlab = expression(alpha),
         ylab = "logLik", ...)
  opar <- par("mar")
  par(mar = c(5, 5, 2, 2))
  matplot(object$alphas, object$cfx, type = "l", xlab = expression(alpha),
          ylab = expression(hat(beta)[j](alpha)), ...)
  text(x = min(object$alphas) + 0.05 * abs(diff(range(object$alphas))),
       y = object$cfx[which.min(object$alphas), ],
       labels = colnames(object$cfx), ...)
  par(mar = opar)
}

.check_bounds <- function(mi, ma) {
  if (ma > mi && ma >= 0 && ma <= 1 && mi <= 1 && mi >= 0)
    ret <- TRUE
  else
    ret <- FALSE
  return(ret)
}
