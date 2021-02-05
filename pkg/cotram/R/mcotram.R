
### mmlt function for count case
mcotram <- function(..., formula = ~ 1, data, theta = NULL,
                    control.outer = list(trace = FALSE), scale = FALSE) {
  
  call <- match.call()
  
  m <- list(...)
  J <- length(m)
  
  ### weights are not yet allowed
  w <- unique(do.call("c", lapply(m, weights)))
  stopifnot(isTRUE(all.equal(w, 1)))
  
  ### warning for todistr != "normal"
  if (any(sapply(m, function(x) x$todistr$name != "normal")))
    warning("One of the models has a non-normal inverse link function F_Z. ML
              optimization still works but no interpretation in the
              Gaussian copula framework is possible, though the lambdas still serve
              as coefficients for the transformation functions.")
  ### check if data is count vector
  lu <- lapply(m, function(mod) {
    iY <- get("iY", environment(mod$parm))
    if (is.null(iY)) 
      stop("not a count outcome.")
    fixed <- get("fixed", environment(mod$parm))
    offset <- get("offset", environment(mod$parm))
    tmp <- attr(iY$Yleft, "constraint")
    wf <- !colnames(iY$Yleft) %in% names(fixed)
    iY$Yleft <- iY$Yleft[, wf,drop = FALSE]
    iY$Yright <- iY$Yright[, wf,drop = FALSE]
    attr(iY$Yleft, "constraint") <- tmp
    attr(iY$Yright, "constraint") <- tmp
    list(lower = iY$Yleft, upper = iY$Yright)
  })
  
  Jp <- J * (J - 1) / 2
  
  bx <- formula
  if (inherits(formula, "formula"))
    bx <- as.basis(formula, data)
  lX <- model.matrix(bx, data = data)
  
  N <- nrow(lX)
  nobs <- sapply(lu, function(m) nrow(m$lower))
  stopifnot(length(unique(nobs)) == 1L)
  
  npar <- sum(sapply(lu, function(m) ncol(m$lower))) + Jp * ncol(lX)
  
  Ylower <- do.call("bdiag", lapply(lu, function(m) m$lower))
  Yupper <- do.call("bdiag", lapply(lu, function(m) m$upper))
  
  cnstr <- do.call("bdiag", lapply(lu, function(m) attr(m$lower, "constraint")$ui))
  ui <- bdiag(cnstr, Diagonal(Jp * ncol(lX)))
  ci <- do.call("c", lapply(lu, function(m) attr(m$lower, "constraint")$ci))
  ci <- c(ci, rep(-Inf, Jp * ncol(lX)))
  ui <- ui[is.finite(ci),]
  ci <- ci[is.finite(ci)]
  ui <- as(ui, "matrix")
  
  idx <- 1
  S <- 1
  if (J > 2) {
    S <- matrix(rep(rep(1:0, (J - 1)), c(rbind(1:(J - 1), Jp))), nrow = Jp)[, -J]
    idx <- unlist(lapply(colSums(S), seq_len))
  }
  
  ### catch constraint violations here
  .log <- function(x) {
    return(log(pmax(.Machine$double.eps, x)))
  }
  
  
  ll <- function(par) {
    
    mpar <- par[1:ncol(Ylower)]
    cpar <- matrix(par[-(1:ncol(Ylower))], nrow = ncol(lX))
    
    ### Ylower = NaN means -Inf and Yupper == NaN means +Inf
    ### (corresponding to probabilities 0 and 1)
    Yp_l <- matrix(Ylower %*% mpar, nrow = N)
    Yp_l[is.na(Yp_l)] <- -Inf
    Yp_u <- matrix(Yupper %*% mpar, nrow = N)
    Yp_u[is.na(Yp_u)] <- Inf
    
    Xp <- lX %*% cpar
    
    A <- Yp_u[, idx] * Xp
    B_l <- A %*% S + Yp_l[,-1]
    B_u <- A %*% S + Yp_u[,-1]
    C_l <- cbind(Yp_l[,1], B_l)
    C_u <- cbind(Yp_u[,1], B_u)
    
    ret <- 0
    for (j in 1:J) {
      F_Zj <- m[[j]]$todistr$p
      ret <- ret + sum(.log(F_Zj(C_u[, j]) - F_Zj(C_l[, j])))
    }
    
    return(-ret)
  }
  
  sc <- function(par) {
    
    mpar <- par[1:ncol(Ylower)]
    cpar <- matrix(par[-(1:ncol(Ylower))], nrow = ncol(lX))
   
    ### Ylower = NaN means -Inf and Yupper == NaN means +Inf
    ### (corresponding to probabilities 0 and 1)
    Yp_l <- matrix(Ylower %*% mpar, nrow = N)
    Yp_l[is.na(Yp_l)] <- -Inf
    Yp_u <- matrix(Yupper %*% mpar, nrow = N)
    Yp_u[is.na(Yp_u)] <- Inf
    Xp <- lX %*% cpar
    
    A <- Yp_u[, idx] * Xp
    B_l <- A %*% S + Yp_l[,-1]
    B_u <- A %*% S + Yp_u[,-1]
    C_l <- cbind(Yp_l[,1], B_l)
    C_u <- cbind(Yp_u[,1], B_u)
    
    C1 <- C_l
    for (j in 1:J) {
      f_Zj <- m[[j]]$todistr$d
      F_Zj <- m[[j]]$todistr$p
      C1[, j] <- (f_Zj(C_u[, j]) - f_Zj(C_l[, j]))/(F_Zj(C_u[, j]) - F_Zj(C_l[, j]))
    }
    
    L <- diag(0, J)
    L[upper.tri(L)] <- 1:Jp
    L <- t(L)
    
    mret <- vector(length = J, mode = "list")
    for (k in 1:J) {
      f_Zk <- m[[k]]$todistr$d
      F_Zk <- m[[k]]$todistr$p
      lu[[k]]$lower[is.infinite(lu[[k]]$lower)] <- 0
      mret[[k]] <- colSums((f_Zk(C_u[, k])*lu[[k]]$upper - f_Zk(C_l[, k])*lu[[k]]$lower)/
                             (F_Zk(C_u[, k]) - F_Zk(C_l[, k])))
      Lk <- L[,k]
      D <- cbind(matrix(rep(0, k*N), nrow = N), Xp[,Lk[Lk > 0]])
      mret[[k]] <- mret[[k]] + colSums(rowSums(C1 * D) * lu[[k]]$upper)
    }
    
    cret <- vector(length = J - 1, mode = "list")
    for (k in 1:(J - 1)) { # go over rows
      C2 <- matrix(rep(C1[, k+1], k), ncol = k)
      tmp <- C2 * Yp_u[,1:k]
      ret <- c()
      l <- ncol(lX)
      for (i in 1:k) {
        tmp1 <- matrix(rep(tmp[,i], l), ncol = l)
        ret <- c(ret, colSums(tmp1 * lX))
      }
      cret[[k]] <- ret
    }
    
    mret <- -do.call("c", mret)
    cret <- -do.call("c", cret)
    c(mret, cret)
  }
  
  ### user-defined starting parameters for optimization
  if(!is.null(theta)) {
    start <- unname(theta)
  }
  else {
    # if(inherits(formula, "formula") && formula == ~1) {
      ### don't bother with .start(), simply use the marginal coefficients
      ### and zero for the lambda parameters
      start <- do.call("c", lapply(m, function(mod) coef(as.mlt(mod))))
      start <- c(start, rep(0, Jp * ncol(lX)))
    # }
    # else { # formula != ~ 1
      # start <- .start(m, bx = bx, data = data)
      # start <- c(start$mpar, c(t(start$cpar)))
    }
  
  
  ### does this work for count too? m$lower or upper?
  if(scale) {
    warning("scale = TRUE not yet implemented. setting scale = FALSE")
    scale <- FALSE
  }
  # if (scale) {
  #   Ytmp <- cbind(do.call("cbind", lapply(lu, function(m) m$lower)), 
  #                 kronecker(matrix(1, ncol = Jp), lX))
  #   Ytmp[!is.finite(Ytmp)] <- NA
  #   scl <- apply(abs(Ytmp), 2, max, na.rm = TRUE)
  #   lt1 <- scl < 1.1
  #   gt1 <- scl >= 1.1
  #   scl[gt1] <- 1 / scl[gt1]
  #   scl[lt1] <- 1
  #   start <- start / scl
  #   if (!is.null(ui))
  #     ui <- t(t(ui) * scl)
  #   f <- function(par) ll(scl * par, ui = ui, ci = ci)
  #   g <- function(par) sc(scl * par) * scl
  # } else {
    f <- function(par) ll(par)
    g <- sc
  # }

  opt <- alabama::auglag(par = start, fn = f, gr = g,
                         hin = function(par) ui %*% par - ci, 
                         hin.jac = function(par) ui,
                         control.outer = control.outer)[c("par", 
                                                          "value", 
                                                          "gradient",
                                                          "hessian")]
  if (scale) opt$par <- opt$par * scl
  
  opt$ll <- ll
  opt$sc <- sc
  opt
  mpar <- opt$par[1:(sum(sapply(lu, function(m) ncol(m$lower))))]

  mlist <- split(mpar, sf <- rep(factor(1:J), sapply(lu, function(m) ncol(m$lower))))
  mmod <- vector(mode = "list", length = J)
  for (j in 1:J) {
    mmod[[j]] <- as.mlt(m[[j]])
    coef(mmod[[j]]) <- mlist[[j]]
  }
  cpar <- matrix(opt$par[-(1:length(mpar))], ncol = Jp)
  
  gaussian <- all.equal("normal", unique(sapply(mmod, function(x) x$todistr$name)))
  
  nm <- abbreviate(sapply(m, function(x) x$model$response), 4)
  lnm <- matrix(paste0(matrix(nm, nrow = J, ncol = J), ".",
                       matrix(nm, nrow = J, ncol = J, byrow = TRUE)), nrow = J)
  cnm <- paste0(rep(lnm[lower.tri(lnm)], each = ncol(lX)), ".", rep(colnames(lX), Jp))
  names(opt$par) <- c(paste0(nm[sf], ".", do.call("c", lapply(mlist, names))), cnm)
  
  ret <- list(marginals = mmod, formula = formula, bx = bx, data = data,
              call = call,
              gaussian = gaussian,
              pars = list(mpar = mpar, cpar = cpar),
              par = opt$par, ll = ll, sc = sc, logLik = opt$value,
              hessian = opt$hessian)
  class(ret) <- c("mcotram", "mmlt")
  ret
}

predict.mcotram <- function(object, newdata = object$data, marginal = 1L,
                            type = c("trafo", "distribution", "density"), ...) {
  type <- match.arg(type)
  if (!object$gaussian & marginal != 1L)
    stop("Cannot compute marginal distribution from non-gaussian joint model")
  
  ### predicting marginal transformation functions
  ret <- lapply(object$marginals[marginal], function(m)
    predict.cotram(m, newdata = newdata, type = "trafo", ...))
  Vx <- coef(object, newdata = newdata, type = "Sigma")

  if (type == "distribution") {
    ret <- lapply(1:length(ret), function(i) {
      tmp <- t(t(ret[[i]]) / sqrt(Vx$diag[,marginal]))
      pnorm(tmp)
    })
  }
  if (type == "density") {
    newdata_m1 <- newdata
    y <- unlist(lapply(object$marginals[marginal], function(m)
      variable.names(m, "response")))
    if (y %in% names(newdata_m1)) newdata_m1[,y] <- newdata_m1[,y] - 1L
    ret_m1 <- lapply(object$marginals[marginal], function(m)
      predict(m, newdata = newdata_m1, ...))
    ret <- lapply(1:length(ret), function(i) {
      tmp <- t(t(ret[[i]]) / sqrt(Vx$diag[,marginal]))
      tmp_m1 <- t(t(ret_m1[[i]]) / sqrt(Vx$diag[,marginal]))
      tmp_m1[is.na(tmp_m1)] <- -Inf
      pnorm(tmp) - pnorm(tmp_m1)
    })
  }
  if (length(ret) == 1) return(ret[[1]])
  ret
}


# coef.mmlt <- function(object, newdata = object$data, 
#                       type = c("all", "marginal", "Lambda", "Lambdainv", "Sigma", "Corr"), 
#                       ...)
# {
#   
#   type <- match.arg(type)
#   if (type == "all") return(object$par)
#   if (type == "marginal") return(lapply(object$marginals, coef))
#   
#   X <- model.matrix(object$bx, data = newdata)
#   ret <- X %*% object$pars$cpar
#   
#   if (!object$gaussian & type != "Lambda")
#     warning("return value of Lambda() has no direct interpretation")
#   
#   return(switch(type, "Lambda" = ret,
#                 "Lambdainv" = .Solve2(ret),
#                 "Sigma" = .Crossp(.Solve2(ret)),
#                 "Corr" = {
#                   ret <- .Crossp(.Solve2(ret))
#                   isd <- sqrt(ret$diagonal)
#                   if (!is.matrix(isd)) isd <- matrix(isd, nrow = 1)
#                   SS <- c()
#                   J <- length(object$marginals)
#                   for (j in 1:J)
#                     SS <- cbind(SS, isd[,j] * isd[,-(1:j), drop = FALSE])
#                   ret$lower / SS
#                 }))
# }
# 
# summary.mmlt <- function(object, ...) {
#   ret <- list(call = object$call,
#               #                tram = object$tram,
#               test = cftest(object, parm = names(coef(object, with_baseline = FALSE))),
#               ll = logLik(object))
#   class(ret) <- "summary.mmlt"
#   ret
# }
# 
# print.summary.mmlt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
#   cat("\n", "Multivariate conditional transformation model", "\n")
#   cat("\nCall:\n")
#   print(x$call)
#   cat("\nCoefficients:\n")
#   pq <- x$test$test
#   mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
#   colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#   sig <- .Machine$double.eps
#   printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
#                P.values = TRUE, eps.Pvalue = sig)
#   cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), ")", sep = "")
#   cat("\n\n")
#   invisible(x)
# }
# 
# print.mmlt <- function(x, ...) {
#   cat("\n", "Multivariate count conditional transformation model", "\n")
#   cat("\nCall:\n")
#   print(x$call)
#   cat("\nCoefficients:\n")
#   print(coef(x))
#   invisible(x)
# }