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


mmlt <- function(..., formula = ~ 1, data, theta = NULL, diag = FALSE,
                 control.outer = list(trace = FALSE), scale = FALSE,
                 tol = sqrt(.Machine$double.eps)) {
  
  call <- match.call()
  # stopifnot(diag)
  
  m <- lapply(list(...), function(x) as.mlt(x))
  J <- length(m)
  
  ### weights are not yet allowed
  w <- unique(do.call("c", lapply(m, weights)))
  stopifnot(isTRUE(all.equal(w, 1)))
  
  ### warning for todistr != "normal"
  if (any(sapply(m, function(x) x$todistr$name != "normal")))
    warning("One of the models has a non-normal inverse link function F_Z. ML
              optimization still works but has been implemented differently than
              described in the MCTM paper. Hence, no interpretation in the
              Gaussian copula framework is possible, though the lambdas still serve
              as coefficients for the transformation functions.")
  ### check if data is continuous and branch to discrete version here
  
  lu <- lapply(m, function(mod) {
    eY <- get("eY", environment(mod$parm))
    if (is.null(eY)) 
      stop("only continuous outcomes implemented so far.")
    fixed <- get("fixed", environment(mod$parm))
    offset <- get("offset", environment(mod$parm))
    tmp <- attr(eY$Y, "constraint")
    wf <- !colnames(eY$Y) %in% names(fixed)
    eY$Y <- eY$Y[, wf,drop = FALSE]
    attr(eY$Y, "constraint") <- tmp
    list(exact = eY$Y, prime = eY$Yprime)
  })
  
  Jp <- J * (J - 1) / 2 + diag * J
  
  bx <- formula
  if (inherits(formula, "formula"))
    bx <- as.basis(formula, data)
  lX <- model.matrix(bx, data = data)
  
  N <- nrow(lX)
  nobs <- sapply(lu, function(m) nrow(m$exact))
  stopifnot(length(unique(nobs)) == 1L)
  
  Y <- do.call("bdiag", lapply(lu, function(m) m$exact))
  Yprime <- do.call("bdiag", lapply(lu, function(m) m$prime))
  
  cnstr <- do.call("bdiag", 
      lapply(lu, function(m) attr(m$exact, "constraint")$ui))
  ui <- bdiag(cnstr, Diagonal(Jp * ncol(lX)))
  ci <- do.call("c", lapply(lu, function(m) attr(m$exact, "constraint")$ci))

  L <- diag(rep(NA, J))
  L[lower.tri(L, diag = diag)] <- 1:Jp
  di <- diag(L)
  di <- di[!is.na(di)]

  CP <- matrix(1:(Jp*ncol(lX)), nrow = ncol(lX))
  dintercept <- CP[1L, di]
  tci <- rep(-Inf, Jp * ncol(lX))
  tci[dintercept] <- 1 - tol
  D <- Diagonal(Jp * ncol(lX))[dintercept,]
  NL <- Matrix(0, nrow = length(dintercept), ncol = ncol(cnstr))
  ui <- rbind(ui, cbind(NL, -D))

  ci <- c(ci, tci, rep(-1 + tol, length(dintercept)))
  ui <- ui[is.finite(ci),]
  ci <- ci[is.finite(ci)]
  ui <- as(ui, "matrix")
  
  idx <- idx_d <- 1
  S <- 1
  if (J > 2 && !diag) {
    S <- matrix(rep(rep(1:0, (J - 1)), c(rbind(1:(J - 1), Jp))), nrow = Jp)[, -J]
    idx <- unlist(lapply(colSums(S), seq_len))
  }
  
  if (diag) {
    S <- matrix(rep(rep(1:0, J),
                    c(rbind(1:J, Jp))), nrow = Jp)[, -(J + 1)]
    idx <- unlist(lapply(colSums(S), seq_len))
    idx_d <- cumsum(unlist(lapply(colSums(S), sum)))
  }
  
  
  ### catch constraint violations here
  .log <- function(x) {
    return(log(pmax(.Machine$double.eps, x)))
    pos <- (x > .Machine$double.eps)
    if (all(pos)) return(log(x))
    ret[pos] <- log(x[pos])
    return(ret)
  }
  
  if(diag){
    ll <- function(par) {
      
      mpar <- par[1:ncol(Y)]
      cpar <- matrix(par[-(1:ncol(Y))], nrow = ncol(lX))
      
      Yp <- matrix(Y %*% mpar, nrow = N)
      Yprimep <- matrix(Yprime %*% mpar, nrow = N)
      Xp <- lX %*% cpar
      
      A <- Yp[, idx] * Xp
      B <- A %*% S
      
      ret <- sum(.log(Yprimep)) + sum(.log(Xp[, idx_d]))
      for (j in 1:J) {
        ret <- ret + sum(m[[j]]$todistr$d(B[, j], log = TRUE))
      }
      
      return(-ret)
    }
    
    sc <- function(par) {
      
      mpar <- par[1:ncol(Y)] 
      cpar <- matrix(par[-(1:ncol(Y))], nrow = ncol(lX))
      
      Yp <- matrix(Y %*% mpar, nrow = N)
      Yprimep <- matrix(Yprime %*% mpar, nrow = N)
      Xp <- lX %*% cpar
      
      ## this is row-wise
      L <- diag(0, J)
      L[!lower.tri(L)] <- 1:Jp
      L <- t(L)
      
      A <- Yp[, idx] * Xp
      B <- A %*% S
      C <- B
      for (j in 1:J) {
        C[, j] <- m[[j]]$todistr$dd2d(B[, j])
      }
      
      mret <- vector(length = J, mode = "list")
      for (k in 1:J) {
        Lk <- L[,k]
        D <- cbind(matrix(rep(0, (k-1)*N), nrow = N), Xp[,Lk[Lk > 0]])
        mret[[k]] <- colSums(rowSums(C * D) * lu[[k]]$exact) +
          colSums(lu[[k]]$prime / Yprimep[,k])
      }
      
      cret <- vector(length = J, mode = "list")
      for (k in 1:J) {  # go over rows, this is \tilde{k} in the formula
        om_Zk <- m[[k]]$todistr$dd2d
        B1 <- matrix(rep(B[,k], k), ncol = k)
        tmp <- om_Zk(B1) * Yp[,1:k]
        ret <- c()
        l <- ncol(lX)
        if(k > 1) {
          for (i in 1:(k-1)) { ## this is k in formula
          tmp1 <- matrix(rep(tmp[,i], l), ncol = l)
          ret <- c(ret, colSums(tmp1 * lX))
          }
        }
        
        tmp1 <- matrix(rep(tmp[,k], l), ncol = l)
        tmp2 <- matrix(rep(1/Xp[,idx_d[k]], l), ncol = l)
        ret <- c(ret, colSums(tmp1 * lX) + colSums(tmp2 * lX))
        
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
      ### this is never the case because diag happens whenever formula is not ~1
      # if(inherits(formula, "formula") && formula == ~1) {
      #   ### don't bother with .start(), simply use the marginal coefficients
      #   ### and zero for the lambda parameters
      #   start <- do.call("c", lapply(m, function(mod) coef(as.mlt(mod))))
      #   start <- c(start, rep(0, Jp * ncol(lX)))
      # }
      # else { # formula != ~ 1
      start <- .start(m, bx = bx, data = data)
      
      ### FIXME: start$cpar needs to include starting values for diagonal
      ### elements as well!
      # start <- c(start$mpar, c(t(start$cpar)))
      CS <- matrix(0, nrow = ncol(lX), ncol = Jp)
      CS[1L, di] <- 1
      cstart <- c(CS)
      # start <- c(start$mpar, rep(0, Jp*ncol(lX)))
      start <- c(start$mpar, cstart)
      
      # }
    }
  }
  else {
    ll <- function(par) {
      
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
        om_Zk <- m[[k+1]]$todistr$dd2d
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
    if(!is.null(theta)) {
      start <- unname(theta)
    }
    else {
      if(inherits(formula, "formula") && formula == ~1) {
        ### don't bother with .start(), simply use the marginal coefficients
        ### and zero for the lambda parameters
        start <- do.call("c", lapply(m, function(mod) coef(as.mlt(mod))))
        start <- c(start, rep(0, Jp * ncol(lX)))
      }
      else { # formula != ~ 1
        start <- .start(m, bx = bx, data = data)
        start <- c(start$mpar, c(t(start$cpar)))
      }
    }
  }
  
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
    f <- function(par) ll(scl * par)
    g <- function(par) sc(scl * par) * scl
  } else {
    f <- function(par) ll(par)
    g <- sc
  }
  
  if(diag) {
    opt <- alabama::auglag(par = start, fn = f, gr = g,
                           hin = function(par) ui %*% par - ci, 
                           hin.jac = function(par) ui,
                           control.outer = control.outer)[c("par", 
                                                            "value", 
                                                            "gradient",
                                                            "hessian")]
  } else {
    opt <- alabama::auglag(par = start, fn = f, gr = g,
                         hin = function(par) ui %*% par - ci, 
                         hin.jac = function(par) ui,
                         control.outer = control.outer)[c("par", 
                                                          "value", 
                                                          "gradient",
                                                          "hessian")]
  }
  
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
  
  if(diag) {### adapt! make sure that order is right! I think here we use columns
    lnm <- matrix(paste0(matrix(nm, nrow = J, ncol = J), ".",
                       matrix(nm, nrow = J, ncol = J, byrow = TRUE)), nrow = J)
    cnm <- paste0(rep(t(lnm[!upper.tri(lnm)]), each = ncol(lX)), ".",
                rep(colnames(lX), Jp))
  }
  else {
    lnm <- matrix(paste0(matrix(nm, nrow = J, ncol = J), ".",
                       matrix(nm, nrow = J, ncol = J, byrow = TRUE)), nrow = J)
    cnm <- paste0(rep(lnm[lower.tri(lnm)], each = ncol(lX)), ".", rep(colnames(lX), Jp))
    
  }
  
  ## not sure if I should adapt this too?
  names(opt$par) <- c(paste0(nm[sf], ".", do.call("c", lapply(mlist, names))), cnm)
  
  ret <- list(marginals = mmod, formula = formula, bx = bx, data = data,
              call = call,
              gaussian = gaussian, diag = diag,
              pars = list(mpar = mpar, cpar = cpar),
              par = opt$par, ll = ll, sc = sc, logLik = -opt$value,
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

### solve lower triangular matrix (in vector form, saved column-wise!)
### rowwise applicable to matrices
.Solve2 <- function(x, diag) {
  if (!is.matrix(x)) x <- matrix(x, nrow = 1)
  n <- J <- (1 + sqrt(1 + 4 * 2 * ncol(x))) / 2 - diag
  
  if (diag) {
    Jp <- ncol(x)
    
    L <- diag(0, J)
    L[!upper.tri(L)] <- 1:Jp  ## column-wise
    
    idx_d <- diag(L)
    idx_l <- L[lower.tri(L)] ## column-wise
    
    x_orig <- x
    x_diag <- x[, idx_d]
    x <- x[, -idx_d]
    if (!is.matrix(x_diag)) x_diag <- matrix(x_diag, nrow = 1)
    if (!is.matrix(x)) x <- matrix(x, nrow = 1)
    
    idx_f <- rep(1, J-1)
    if(J > 2) {
      for (j in 2:J) {
        idx_f <- c(idx_f, rep(j, J-j))
      }
    }
    
    x <- x / x_diag[, idx_f]
    if(J == 2) x <- t(x)
  }
  
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
  if (!diag) {
    ret
  } else {
    x_norm <- ret
    
    idx_f1 <- 2:J
    if(J > 2) {
      for (j in 3:J) {
        idx_f1 <- c(idx_f1, j:J)
      }
    }
    
    ret <- x_orig
    ret[, idx_l] <- x_norm / x_diag[, idx_f1]
    ret[, idx_d] <- 1/x_diag
    
    ret
  }
} 

# we will use as input Solve2(Xp)
.Crossp <- function(Linv, diag) {
  
  if(!is.matrix(Linv)) Linv <- matrix(Linv, nrow = 1)
  
  # 1 observation
  N <- nrow(Linv)
  Jp <- ncol(Linv)
  J <- (1 + sqrt(1 + 4 * 2 * Jp)) / 2 - diag
  if (N == 1) {
    if (!diag) {
      L <- diag(J)
      L[upper.tri(L)] <- Linv
      L <- t(L)
    }
    else {
      L <- diag(0, J)
      L[!upper.tri(L)] <- Linv
    }
    
    tcp <- tcrossprod(L)
    S_diag <- diag(tcp)
    S_low <- tcp[lower.tri(tcp)]
  }
  # more than 1 observation
  else{
    if (!diag) {
      # J = 1
      S_diag <- 1
      
      if (J == 2) {
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
    } else { ## diag = TRUE
      # J = 1
      S_diag <- Linv^2
      
      if (J == 2) {
        S_diag <- cbind(Linv[, 1]^2, Linv[, 2]^2 + Linv[, 3]^2)
      }
      S_low <- Linv[, 1]*Linv[, 2]
      ## J > 2
      if (J > 2) {
        L <- diag(0, J)
        L[!upper.tri(L)] <- 1:Jp
        
        S <- matrix(rep(rep(1:0, J), c(rbind(1:J, Jp))), nrow = Jp)[, -(J+1)]
        idx_d <- cumsum(unlist(lapply(colSums(S), sum)))
        
        S_low <- Linv*0 ## ensures right length
        for (i in J:1) { # row
          for (j in i:1) { # column
            for (k in 1:j) { # produkt-summanden
              S_low[, L[i, j]] <- S_low[, L[i, j]] + Linv[, L[i, k]]*Linv[, L[j, k]]
            }
          }
        }
        S_diag <- S_low[, idx_d]
        S_low <- S_low[, -idx_d]
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
  
  diagg <- object$diag

  if (!object$gaussian & type != "Lambda")
    warning("return value of Lambda() has no direct interpretation")
  
  return(switch(type, "Lambda" = ret,
                "Lambdainv" = .Solve2(ret, diagg),
                "Sigma" = {
                  .Crossp(.Solve2(ret, diagg), diagg)
                  },
                "Corr" = {
                  ret <- .Crossp(.Solve2(ret, diagg), diagg)
                  isd <- sqrt(ret$diagonal)
                  if (!is.matrix(isd)) isd <- matrix(isd, nrow = 1)
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

print.mmlt <- function(x, ...) {
  cat("\n", "Multivariate conditional transformation model", "\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x))
  if (x$diag) {
     cat("\nDiagonal:\n", "elements are estimated.\n")
  } else { 
    cat("\nDiagonal:\n", "elements are constrained to 1.\n")
    }
  invisible(x)
}
