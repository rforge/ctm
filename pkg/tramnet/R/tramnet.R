# tramnet main function

tramnet <-
  function(model, x, lambda, alpha, constraints = NULL, ...) {
    ### <FIXME> handle offset and maybe fixed parameters </FIXME>
    .tramnet_checks(model = model, x = x, lambda = lambda, alpha = alpha)
    call <- match.call()
    trdat <- .get_tram_data(model)
    stopifnot(trdat$nobs == nrow(x))
    nth <- trdat$npar
    nb <- ncol(x)
    theta <- Variable(nth)
    beta <- Variable(nb)
    prob <- .tramnet_objective(trdat, x, theta, beta,
                               alpha, lambda, constraints)
    res <- solve(prob, ...)
    ret <- list(call = call, model = model, x = x, result = res,
                beta = res$getValue(beta), theta = res$getValue(theta),
                tuning_parm = c(lambda = lambda, alpha = alpha))
    names(ret$beta) <- colnames(x)
    class(ret) <- .tramnet_class(model)
    return(ret)
  }

.tramnet_checks <- function(model, x, lambda, alpha) {
  if (!(inherits(model, "tram") | inherits(model, "mlt")))
    stop("The provided model should be of class 'tram'")
  if (!inherits(x, "matrix"))
    stop("x should be a matrix")
  if (!(lambda >= 0 && alpha >= 0 && alpha <= 1))
    stop("lambda >= 0 and 0 <= alpha <= 1 is required")
  scx <- scale(x)
  centers <- attr(scx, "scaled:center")
  scales <- attr(scx, "scaled:scale")
  condcen <- any(centers > 1e-6)
  condsca <- any(scales > 1 + 1e-6) | any(scales < 1 - 1e-6)
  if ((condcen | condsca) & (lambda > 0))
    warning("Unscaled design matrices are not sensible under regularization.
            Consider scaling and centering the design matrix.")
}

.tramnet_objective <-
  function(trdat, x, theta, beta, alpha, lambda, constraints) {
    xe <- x[trdat$exact$which, , drop = FALSE]
    xl <- x[trdat$censl$which, , drop = FALSE]
    xr <- x[trdat$censr$which, , drop = FALSE]
    xi <- x[trdat$censi$which, , drop = FALSE]
    neg <- ifelse(trdat$negative, -1, 1)
    switch(
      trdat$errdistr,
      "minimum extreme value" = {
        lle <- if (nrow(trdat$exact$ay) == 0)
          0
        else {
          z <- trdat$exact$ay %*% theta + neg * xe %*% beta
          sum_entries(z - exp(z) + log(trdat$exact$aypr %*% theta))
        }
        lll <- if (nrow(trdat$censl$ay) == 0)
          0
        else {
          z <- trdat$censl$ay %*% theta + neg * xl %*% beta
          ### <FIXME> does not work
          sum_entries(log1p(-exp(-exp(z))))
          ### </FIXME>
        }
        llr <- if (nrow(trdat$censr$ay) == 0)
          0
        else {
          z <- trdat$censr$ay %*% theta + neg * xr %*% beta
          sum_entries(-exp(z))
        }
        lli <- if (nrow(trdat$censi$ayl) == 0)
          0
        else {
          xb <- neg * xi %*% beta
          zl <- trdat$censi$ayl %*% theta + xb
          zr <- trdat$censi$ayr %*% theta + xb
          ### <FIXME> does not work
          sum_entries(log(exp(-exp(zr)) - exp(-exp(zl))))
          ### </FIXME>
        }
        ## TODO: Truncation
      },
      "maximum extreme value" = {
        lle <- if (nrow(trdat$exact$ay) == 0)
          0
        else {
          z <- trdat$exact$ay %*% theta + neg * xe %*% beta
          sum_entries(-exp(-z) - z + log(trdat$exact$aypr %*% theta))
        }
        lll <- if (nrow(trdat$censl$ay) == 0)
          0
        else {
          z <- trdat$censl$ay %*% theta + neg * xl %*% beta
          sum_entries(-exp(-z))
        }
        llr <- if (nrow(trdat$censr$ay) == 0)
          0
        else {
          z <- trdat$censr$ay %*% theta + neg * xr %*% beta
          ### <FIXME> does not work
          sum_entries(log1p(-exp(-exp(-z))))
          ### </FIXME>
        }
        lli <- if (nrow(trdat$censi$ayl) == 0)
          0
        else {
          xb <- neg * xi %*% beta
          zl <- trdat$censi$ayl %*% theta + xb
          zr <- trdat$censi$ayr %*% theta + xb
          ### <FIXME> does not work
          sum_entries(log(exp(-exp(-zr)) - exp(-exp(-zl))))
          ### </FIXME>
        }
        ## TODO: Truncation
      },
      "logistic" = {
        lle <- if (nrow(trdat$exact$ay) == 0)
          0
        else {
          z <- trdat$exact$ay %*% theta + neg * xe %*% beta
          sum_entries(-logistic(-z) - logistic(z) + log(trdat$exact$aypr %*% theta))
        }
        lll <- if (nrow(trdat$censl$ay) == 0)
          0
        else {
          z <- trdat$censl$ay %*% theta + neg * xl %*% beta
          sum_entries(-logistic(-z))
        }
        llr <- if (nrow(trdat$censr$ay) == 0)
          0
        else {
          z <- trdat$censr$ay %*% theta + neg * xr %*% beta
          sum_entries(-logistic(z))
        }
        lli <- if (nrow(trdat$censi$ayl) == 0)
          0
        else {
          xb <- neg * xi %*% beta
          zl <- trdat$censi$ayl %*% theta + xb
          zr <- trdat$censi$ayr %*% theta + xb
          ### <FIXME> Does not work
          sum_entries(log(exp(-logistic(-zr)) - exp(-logistic(-zl))))
          ### </FIXME>
        }
        ## TODO: Truncation
      },
      "normal" = {
        if (nrow(trdat$exact$ay) != trdat$nobs)
          stop("BoxCox tramnet not able to deal with any form of censoring")
        z <- trdat$exact$ay %*% theta + neg * xe %*% beta
        const <- -log(sqrt(2 * pi)) * trdat$nobs
        lle <-
          const - sum(z ^ 2 / 2) + sum_entries(log(trdat$exact$aypr %*% theta))
        lll <- llr <- lli <- 0
      }
    )
    obj <- -(lle + lll + llr + lli) +
      lambda * (0.5 * (1 - alpha) * power(p_norm(beta, 2), 2) +
                  alpha * p_norm(beta, 1))
    const <- list(trdat$const$ui %*% theta >= trdat$const$ci)
    if (!is.null(constraints)) {
      const[[2]] <- constraints[[1]] %*% beta >= constraints[[2]]
    }
    prob <- Problem(Minimize(obj), constraints = const)
    return(prob)
  }


.get_tram_data <- function(mod) {
  ## mod: a tram model, usually y ~ 1, but it can also contain covariates.
  ## In this case their coefficients are not counted in the penalty term.
  ## stopifnot(is.null(mod$shiftcoef))
  dat <- mget(
    c("iY", "eY", "offset"),
    envir = environment(mod$parm),
    ifnotfound = list(NULL)
  )
  out <- list()
  np <- length(mod$par)
  out$npar <- np
  out$nobs <- NROW(dat$iY$Yleft) + NROW(dat$eY$Y)
  ## === Censoring
  if (is.null(dat$iY)) {
    out$censl <- list(ay = matrix(0, ncol = np, nrow = 0), which = NULL)
    out$censr <-
      list(ay = matrix(0, ncol = np, nrow = 0), which = NULL)
    out$censi <- list(
      ayl = matrix(0, ncol = np, nrow = 0),
      ayr = matrix(0, ncol = np, nrow = 0),
      which = NULL
    )
  } else {
    idxr <-
      which(is.finite(dat$iY$Yleft[, 1]) &
              !is.finite(dat$iY$Yright[, 1]))
    idxl <-
      which(!is.finite(dat$iY$Yleft[, 1]) &
              is.finite(dat$iY$Yright[, 1]))
    idxi <-
      which(is.finite(dat$iY$Yleft[, 1]) &
              is.finite(dat$iY$Yright[, 1]))
    out$censl <- list(ay = dat$iY$Yright[idxl, , drop = FALSE],
                      which = dat$iY$which[idxl])
    out$censr <- list(ay = dat$iY$Yleft[idxr, , drop = FALSE],
                      which = dat$iY$which[idxr])
    out$censi <- list(
      ayl = dat$iY$Yleft[idxi, , drop = FALSE],
      ayr = dat$iY$Yright[idxi, , drop = FALSE],
      which = dat$iY$which[idxi]
    )
  }

  ## === Exact observations
  if (is.null(dat$eY)) {
    out$exact <- list(
      ay = matrix(0, ncol = np, nrow = 0),
      aypr = matrix(0, ncol = np, nrow = 0),
      which = NULL
    )
  } else {
    out$exact <-
      list(
        ay = dat$eY$Y,
        aypr = dat$eY$Yprime,
        which = dat$eY$which
      )
  }
  ## === Offsets, weights, error distribution, etc
  out$offset <- dat$offset
  out$weights <- mod$weights
  out$negative <- mod$negative
  out$errdistr <- mod$todistr$name
  ## === Constraints
  if (!is.null(dat$eY)) {
    out$constr <- attr(dat$eY$Y, "constraint")
  } else {
    out$constr <- attr(dat$iY$Yleft, "constraint")
  }
  ## === Truncation
  ## Left
  idxi <- which(is.finite(dat$iY$trunc$left[, 1]) &
                  !is.finite(dat$iY$trunc$right[, 1]))
  idxe <- which(is.finite(dat$eY$trunc$left[, 1]) &
                  !is.finite(dat$eY$trunc$right[, 1]))
  out$truncl <-
    list(
      ay = rbind(dat$iY$trunc$left[idxi, , drop = FALSE],
                 dat$eY$trunc$left[idxe, , drop = FALSE]),
      which = c(dat$iY$which[idxi], dat$eY$which[idxe])
    )
  ## Right
  idxi <- which(!is.finite(dat$iY$trunc$left[, 1]) &
                  is.finite(dat$iY$trunc$right[, 1]))
  idxe <- which(!is.finite(dat$eY$trunc$left[, 1]) &
                  is.finite(dat$eY$trunc$right[, 1]))
  out$truncr <-
    list(
      ay = rbind(dat$iY$trunc$right[idxi, , drop = FALSE],
                 dat$eY$trunc$right[idxe, , drop = FALSE]),
      which = c(dat$iY$which[idxi], dat$eY$which[idxe])
    )
  ## Interval
  idxi <- which(is.finite(dat$iY$trunc$left[, 1]) &
                  is.finite(dat$iY$trunc$right[, 1]))
  idxe <- which(is.finite(dat$eY$trunc$left[, 1]) &
                  is.finite(dat$eY$trunc$right[, 1]))
  out$trunci <-
    list(
      ayl = rbind(dat$iY$trunc$left[idxi, , drop = FALSE],
                  dat$eY$trunc$left[idxe, , drop = FALSE]),
      ayr = rbind(dat$iY$trunc$right[idxi, , drop = FALSE],
                  dat$eY$trunc$right[idxe, , drop = FALSE]),
      which = c(dat$iY$which[idxi], dat$eY$which[idxe])
    )
  return(out)
}

.tramnet_class <- function(model) {
  if (inherits(model, "Lm"))
    return(c("tramnet_Lm", "tramnet"))
  return("tramnet")
}
