### compute more sensible starting values
### this is essentially the old mmlt() function
.start <- function(by, bx = NULL, data, ...) {

    J <- length(by)

    mctm <- vector(mode = "list", length = J)
    mmlt <- vector(mode = "list", length = J)
    mctm[[1]] <- by[[1]]$model
    mmlt[[1]] <- mlt(mctm[[1]], data = data, ...)
    pdat <- data
    htotal <- "~ 1"

    for (j in 2:J) {
        hhat <- paste("hhat", j - 1, sep = "_")
        htotal <- c(htotal, hhat)
        data[[hhat]] <- predict(mmlt[[j - 1]], newdata = pdat, 
                                type = "trafo")
        pdat[[hhat]] <- 0
        bhi <- as.basis(as.formula(paste(htotal, collapse = "+")), 
                        data = data, remove_intercept = TRUE)
        if (!is.null(bx)) {
            shift <- b(bh = bhi, bx = bx)
            if (!is.null(by[[j]]$model$bases$shifting))
                shift <- c(shift = by[[j]]$model$bases$shifting, bhbx = b(bh = bhi, bx = bx))
            mctm[[j]] <- ctm(by[[j]]$model$bases$response, 
                             interacting = by[[j]]$model$bases$interacting,
                             shifting = shift,
                             todistr = "Normal")
        } else {
            shift <- bhi
            if (!is.null(by[[j]]$model$bases$shifting))
                shift <- c(shift = by[[j]]$model$bases$shifting, bhbx = b(bh = bhi, bx = bx))
            mctm[[j]] <- ctm(by[[j]]$model$bases$response, 
                             interacting = by[[j]]$model$bases$interacting,
                             shifting = shift,
                             todistr = "Normal")
        }
        ### set todistr
        mctm[[j]]$todistr <- by[[j]]$todistr
        ### get marginal parameters as starting values
        theta <- coef(mctm[[j]])
        theta[] <- 0
        theta[names(coef(by[[j]]))] <- coef(by[[j]])
        mmlt[[j]] <- mlt(mctm[[j]], data = data, theta = theta, ...)
    }

    ### postprocess parameters
    p <- ncol(model.matrix(bx, data = data))
    cf <- lapply(mmlt, coef)
    mpar <- c()
    for (i in 1:length(cf))
        mpar <- c(mpar, cf[[i]][names(coef(by[[i]]))])
    cpar <- c()
    j <- 1
    for (i in 2:length(cf)) {
        cp <- cf[[i]][grep("hhat", names(cf[[i]]))]
        cpar <- rbind(cpar, matrix(cp, ncol = p))
    }
    list(mpar = mpar, cpar = cpar)
}

# omegas in dd2d argument
mmlt <- function(..., formula = ~ 1, data, theta = NULL,
                 control.outer = list(trace = FALSE), scale = FALSE) {

  call <- match.call()
  
  m <- list(...)
  J <- length(m)

  ### weights are not yet allowed
  w <- unique(do.call("c", lapply(m, weights)))
  stopifnot(isTRUE(all.equal(w, 1)))

  ### check if data is continuous and branch to discrete version here

  lu <- lapply(m, function(mod) {
    eY <- get("eY", environment(mod$parm))
    fixed <- get("fixed", environment(mod$parm))
    offset <- get("offset", environment(mod$parm))
    tmp <- attr(eY$Y, "constraint")
    wf <- !colnames(eY$Y) %in% names(fixed)
    eY$Y <- eY$Y[, wf,drop = FALSE]
    attr(eY$Y, "constraint") <- tmp
    list(exact = eY$Y, prime = eY$Yprime)
  })
  
  Jp <- J * (J - 1) / 2

  bx <- formula
  if (inherits(formula, "formula"))
    bx <- as.basis(formula, data)
  lX <- model.matrix(bx, data = data)

  N <- nrow(lX)
  nobs <- sapply(lu, function(m) nrow(m$exact))
  stopifnot(length(unique(nobs)) == 1L)

  npar <- sum(sapply(lu, function(m) ncol(m$exact))) + Jp * ncol(lX)
  
  Y <- do.call("bdiag", lapply(lu, function(m) m$exact))
  Yprime <- do.call("bdiag", lapply(lu, function(m) m$prime))
  
  cnstr <- do.call("bdiag", lapply(lu, function(m) attr(m$exact, "constraint")$ui))
  ui <- bdiag(cnstr, Diagonal((Jp * ncol(lX))))
  ci <- do.call("c", lapply(lu, function(m) attr(m$exact, "constraint")$ci))
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
    pos <- (x > .Machine$double.eps)
    if (all(pos)) return(log(x))
    ret[pos] <- log(x[pos])
    return(ret)
  }


  ll <- function(par, ui, ci) {

    mpar <- par[1:ncol(Y)]
    cpar <- matrix(par[-(1:ncol(Y))], nrow = ncol(lX))
    
    Yp <- matrix(Y %*% mpar, nrow = N)
    Yprimep <- matrix(Yprime %*% mpar, nrow = N)
    Xp <- lX %*% cpar
    
    A <- Yp[, idx] * Xp
    B <- A %*% S + Yp[,-1]
    C <- cbind(Yp[,1], B)
    
    ret <- sum(.log(Yprimep))
    for (j in 1:J) {
      ret <- ret + sum(m[[j]]$todistr$d(C[, j], log = TRUE))
    }
    
    return(-ret)
  }
  
  sc <- function(par) {
    
    mpar <- par[1:ncol(Y)] 
    cpar <- matrix(par[-(1:ncol(Y))], nrow = ncol(lX))
    
    Yp <- matrix(Y %*% mpar, nrow = N)
    Yprimep <- matrix(Yprime %*% mpar, nrow = N)
    Xp <- lX %*% cpar
    
    L <- diag(0, J)
    L[upper.tri(L)] <- 1:Jp
    L <- t(L)
    
    A <- Yp[, idx] * Xp
    B <- A %*% S + Yp[,-1]
    C <- cbind(Yp[,1], B)
    C1 <- C
    for (j in 1:J) {
      C1[, j] <- m[[j]]$todistr$dd2d(C[, j])
    }
    
    mret <- vector(length = J, mode = "list")
    for (k in 1:J) {
      Lk <- L[,k]
      D <- cbind(matrix(rep(0, (k-1)*N), nrow = N), 1, Xp[,Lk[Lk > 0]])
      mret[[k]] <- colSums(rowSums(C1 * D) * lu[[k]]$exact) +
        colSums(lu[[k]]$prime / Yprimep[,k])
    }
    
    cret <- vector(length = J - 1, mode = "list")
    for (k in 1:(J - 1)) {  # go over rows
      om_Zk <- m[[k]]$todistr$dd2d
      B1 <- matrix(rep(B[,k], k), ncol = k)
      tmp <- om_Zk(B1) * Yp[,1:k]
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
  ### useful for bootstrapping
  if(!is.null(theta)) {
    start <- theta
  }
  else {
    start <- .start(m, bx = bx, data = data)
    start <- c(start$mpar, c(t(start$cpar)))
  }

### this should give the same likelihood as logLik(mmlt()) of the "old"
### version
# print(ll(start))

  if (scale) {
    Ytmp <- cbind(do.call("cbind", lapply(lu, function(m) m$exact)), 
                  kronecker(matrix(1, ncol = Jp), lX))
    Ytmp[!is.finite(Ytmp)] <- NA
    scl <- apply(abs(Ytmp), 2, max, na.rm = TRUE)
    lt1 <- scl < 1.1
    gt1 <- scl >= 1.1
    scl[gt1] <- 1 / scl[gt1]
    scl[lt1] <- 1
    start <- start / scl
    if (!is.null(ui))
        ui <- t(t(ui) * scl)
    f <- function(gamma) ll(scl * gamma, ui = ui, ci = ci)
    g <- function(gamma) sc(scl * gamma) * scl
  } else {
    f <- function(par) ll(par, ui = ui, ci = ci)
    g <- sc
  }

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
  mpar <- opt$par[1:(sum(sapply(lu, function(m) ncol(m$exact))))]

  mlist <- split(mpar, sf <- rep(factor(1:J), sapply(lu, function(m) ncol(m$exact))))
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
  class(ret) <- "mmlt"
  ret
}

predict.mmlt <- function(object, newdata, marginal = 1L, 
                         type = c("trafo", "distribution", "density"), ...) {
    type <- match.arg(type)
    if (!object$gaussian & marginal != 1L)
        stop("Cannot compute marginal distribution from non-gaussian joint model")
    ret <- lapply(object$marginals[marginal], function(m)
        predict(m, newdata = newdata, ...))
    Vx <- coef(object, newdata = newdata, type = "Sigma")
    ### first formula in Section 2.4
    if (type == "distribution") {
        ret <- lapply(1:length(ret), function(i) {
                      tmp <- t(t(ret[[i]]) / sqrt(Vx$diag[,marginal]))
                      pnorm(tmp)
        })
    }
    if (type == "density") {
        hprime <- lapply(object$marginals[marginal], function(m) {
            dr <- 1
            names(dr) <- m$model$response
            predict(m, newdata = newdata, deriv = dr, ...)
        })
        ret <- lapply(1:length(ret), function(i) {
                      tmp <- t(t(ret[[i]]) / sqrt(Vx$diag[,marginal]))
                      t(t(dnorm(tmp)) / sqrt(Vx$diag[,marginal]))  * hprime[[i]]
        })
    }
    if (length(ret) == 1) return(ret[[1]])
    ret
}

### solve lower triangular matrix (in vector form)
### rowwise applicable to matrices
.Solve2 <- function(x) {
    if (!is.matrix(x)) x <- matrix(x, nrow = 1)
    n <- (1 + sqrt(1 + 4 * 2 * ncol(x))) / 2
    xij <- function(x = NULL, i, j) {
        if (i == j) return(1)
        if (j == 1) {
            ret <- i - 1
        } else {
            idx <- n - (1:(n - 1))
            ret <- sum(idx[1:(j - 1)]) + (i - (n - idx[j]))
        }
        if (is.null(x))
            return(ret)
        return(x[,ret])
    }
    ret <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for (i in 2:n) {
        for (j in 1:(i - 1)) {
            s <- 0
            for (k in j:(i - 1))
                s <- s + xij(x, i, k) * xij(ret, k, j)
             ret[, xij(NULL, i, j)] <- -s
        }
    }
    ret
} 

# we will use as input Solve2(Xp)
.Crossp <- function(Linv) {
  # 1 observation
  N <- nrow(Linv)
  Jp <- ncol(Linv)
  J <- (1 + sqrt(1 + 4 * 2 * Jp)) / 2
  if (N == 1) {
    L <- diag(1, J)
    L[upper.tri(L)] <- Linv
    L <- t(L)
    tcp <- tcrossprod(L)
    S_diag <- diag(tcp)
    S_low <- tcp[lower.tri(tcp)]
  }
  # more than 1 observation
  else{
    # J = 1
    S <- S_diag <- 1
    # J = 2
    if(ncol(Linv) == 1) {
      S_diag <- cbind(rep(1, N), matrix(1, nrow = N, ncol = 1) + Linv^2)
      }
    S_low <- Linv

    if (J > 2) {
      L <- diag(0, J)
      L[upper.tri(L)] <- 1:Jp
      L <- t(L)

      S <- matrix(rep(rep(1:0, (J - 1)), c(rbind(1:(J - 1), Jp))), nrow = Jp)[, -J]
      S_diag <- cbind(rep(1, N), matrix(1, nrow = N, ncol = J-1) + Linv^2 %*% S)

      for (i in J:3) { #zeile
        for (j in (i-1):2) { #spalte
          for (k in 1:(j-1)) { #produkt-summanden
            S_low[, L[i, j]] <- S_low[, L[i, j]] + Linv[, L[i, k]]*Linv[, L[j, k]]
          }
        }
      }
    }
  }
  ret <- list(lower = S_low, diagonal = S_diag)
  return(ret)
}

logLik.mmlt <- function(object, ...) {
    ret <- object$logLik
    attr(ret, "df") <- length(object$par)
    class(ret) <- "logLik"
    ret
}

coef.mmlt <- function(object, newdata = object$data, 
                      type = c("all", "marginal", "Lambda", "Lambdainv", "Sigma", "Corr"), 
                      ...)
{

    type <- match.arg(type)
    if (type == "all") return(object$par)
    if (type == "marginal") return(lapply(object$marginals, coef))

    X <- model.matrix(object$bx, data = newdata)
    ret <- X %*% object$pars$cpar

    if (!object$gaussian & type != "Lambda")
        warning("return value of Lambda() has no direct interpretation")

    return(switch(type, "Lambda" = ret,
                        "Lambdainv" = .Solve2(ret),
                        "Sigma" = .Crossp(.Solve2(ret)),
                        "Corr" = {
      ret <- .Crossp(.Solve2(ret))
      isd <- sqrt(ret$diagonal)
      SS <- c()
      J <- length(object$marginals)
      for (j in 1:J)
          SS <- cbind(SS, isd[,j] * isd[,-(1:j), drop = FALSE])
      ret$lower / SS
    }))
}

vcov.mmlt <- function(object, ...) {
    ret <- solve(object$hessian)
    rownames(ret) <- colnames(ret) <- names(coef(object))
    ret
}

summary.mmlt <- function(object, ...) {
    ret <- list(call = object$call,
#                tram = object$tram,
                test = cftest(object, parm = names(coef(object, with_baseline = FALSE))),
                ll = logLik(object))
    class(ret) <- "summary.mmlt"
    ret
}

print.summary.mmlt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\n", "Multivariate conditional transformation model", "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    pq <- x$test$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    sig <- .Machine$double.eps
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
        P.values = TRUE, eps.Pvalue = sig)
    cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), ")", sep = "")
    cat("\n\n")
    invisible(x)
}
