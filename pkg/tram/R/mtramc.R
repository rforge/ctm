
### marginally interpretable linear transformation models for clustered
### observations

mtramc <- function(object, formula, data, standardise = FALSE,
                   grd = SparseGrid::createSparseGrid(type = "KPU", dimension = length(rt$cnms[[1]]), 
                                                      k = 10),
                   Hessian = FALSE,
                   ...) {

    stopifnot(inherits(object, "mlt_fit"))

    bar.f <- lme4::findbars(formula)
    mf <- model.frame(lme4::subbars(formula), data = data)
    rt <- lme4::mkReTrms(bar.f, mf)

    ZtW <- rt$Zt
    Lambdat <- rt$Lambdat
    Lind <- rt$Lind
    mapping <- function(theta)
        theta[Lind]
    theta <- rt$theta

    eY <- get("eY", environment(object$loglik))
    iY <- get("iY", environment(object$loglik))
    fixed <- get("fixed", environment(object$loglik))
    offset <- get("offset", environment(object$loglik))
    wf <- 1:length(coef(as.mlt(object)))
    if (!is.null(eY)) {
        tmp <- attr(eY$Y, "constraint")
        wf <- !colnames(eY$Y) %in% names(fixed)
        eY$Y <- eY$Y[, wf,drop = FALSE]
        attr(eY$Y, "constraint") <- tmp
    }
    if (!is.null(iY)) {
        tmp <- attr(iY$Yleft, "constraint")
        wf <- !colnames(iY$Yleft) %in% names(fixed)
        iY$Yleft <- iY$Yleft[, wf,drop = FALSE]
        iY$Yright <- iY$Yright[, wf,drop = FALSE]
        attr(iY$Yleft, "constraint") <- tmp
        attr(iY$Yright, "constraint") <- tmp
    }
    if (length(eY$which) > 0 && length(iY$which))
        stop("cannot deal with mixed censoring")

    w <- object$weights

    NORMAL <- FALSE
    if (object$todistr$name == "normal") {
        NORMAL <- TRUE
        PF <- function(z) z
    } else {
        P <- object$todistr$p
        PF <- function(z) qnorm(P(z))
    }

    gr <- NULL

    if (length(eY$which) > 0) {
        if (standardise) stop("standardise not yet implemented for continuous case")
        L <- Cholesky(crossprod(Lambdat %*% ZtW), LDL = FALSE, Imult=1)
        ll <- function(parm) {
            theta <- parm[1:ncol(eY$Y)]
            gamma <- parm[-(1:ncol(eY$Y))]
            z <- PF(c(eY$Y %*% theta + offset))
            Lambdat@x[] <- mapping(gamma)
            L <- update(L, t(Lambdat %*% ZtW), mult = 1)
            Linv <- solve(as(L, "Matrix"))
            logdet <- 2 * determinant(L, logarithm = TRUE)$modulus
            ret <- -0.5 * (logdet + sum((Linv %*% z)^2))
            ret <- ret - object$loglik(theta, weights = w) + .5 * sum(z^2)
            return(-ret)
        }
        X <- eY$Y
        if (NORMAL) {
            gr <- function(parm) {
                theta <- parm[1:ncol(eY$Y)]
                gamma <- parm[-(1:ncol(eY$Y))]
                devLambda <- devSigma <- vector(mode = "list", 
                                                length = length(gamma))
                dgamma <- numeric(length(gamma))
                z <- c(eY$Y %*% theta + offset)
                Lambdat@x[] <- mapping(gamma)
                L <- update(L, t(Lambdat %*% ZtW), mult = 1)
                Linv <- solve(as(L, "Matrix"))
                Sigma <- tcrossprod(as(L, "Matrix"))
                SigmaInv <- crossprod(Linv)
                LambdaInd <- t(Lambdat)
                LambdaInd@x[] <- 1:length(gamma)
                for (i in 1:length(gamma)) {
                    ### Wang & Merkle (2018, JSS) compute derivative of G!
                    ### We need derivative of Lambda!
                    dLtL <- (LambdaInd == i) %*% Lambdat
                    devLambda[[i]] <- dLtL + t(dLtL)
                    devSigma[[i]] <- crossprod(ZtW, devLambda[[i]] %*% ZtW)
                    t1 <- SigmaInv %*% devSigma[[i]]
                    t2 <- t1 %*% SigmaInv
                    dgamma[i] <- -.5 * (sum(diag(t1)) - crossprod(z, t2) %*% z)
                }
                dtheta <- - z %*% SigmaInv %*% eY$Y + 
                            colSums(eY$Yprime / c(eY$Yprime %*% theta))
                return(-c(c(as(dtheta, "matrix")), dgamma))
            }
        }
    } else {
        stopifnot(length(rt$flist) == 1)
        grp <- rt$flist[[1]]
        idx <- split(1:length(grp), grp)
        wh <- 1:length(rt$cnms[[1]])
        ### .Marsaglia_1963 expects t(nodes) !!!
        grd$nodes <- t(qnorm(grd$nodes))

## don't spend time on Matrix dispatch
        mZtW <- as(ZtW, "matrix")
        zt <- lapply(idx, function(i) {
            z <- mZtW[,i,drop = FALSE]
            t(z[base::rowSums(abs(z)) > 0,,drop = FALSE])
        })
        mLt <- t(as(Lambdat[wh, wh], "matrix"))
        ONE <- matrix(1, nrow = NCOL(mLt))

## because this needs a lot of time in Matrix
#                z <- ZtW[,i,drop = FALSE]
#                z <- z[rowSums(abs(z)) > 0,,drop = FALSE]
#                V <- t(as(Lambdat[wh, wh] %*% z, "matrix"))

        ll <- function(parm) {
            theta <- parm[1:ncol(iY$Yleft)]
            gamma <- parm[-(1:ncol(iY$Yleft))]
            Lambdat@x[] <- mapping(gamma)
            lplower <- c(iY$Yleft %*% theta + offset)
            lplower[is.na(lplower)] <- -Inf
            lpupper <- c(iY$Yright %*% theta + offset)
            lpupper[is.na(lpupper)] <- Inf

            ret <- sapply(1:length(idx), function(i) {
                V <- zt[[i]] %*% mLt
                ii <- idx[[i]]
                if (standardise) {
                    sd <- c(sqrt((V^2) %*% ONE + 1))
                    zlower <- PF(lplower[ii] / sd) * sd
                    zupper <- PF(lpupper[ii] / sd) * sd
                } else {
                    zlower <- PF(lplower[ii])
                    zupper <- PF(lpupper[ii])
                }
                .Marsaglia_1963(zlower, zupper, mean = 0, V = V, 
                                do_qnorm = FALSE, grd = grd)
            })
            return(-sum(log(pmax(.Machine$double.eps, ret))))
        }
        X <- iY$Yleft
    }            

    ui <- attr(X, "constraint")$ui[, wf, drop = FALSE]
    ci <- attr(X, "constraint")$ci
    ui <- as(bdiag(ui, Diagonal(length(theta))), "matrix")
    ci <- c(ci, rt$lower)

    start <- c(coef(as.mlt(object), fixed = FALSE), theta)

    if (is.null(gr)) {
        opt <- alabama::auglag(par = start, fn = ll, 
            hin = function(par) ui %*% par - ci, 
            hin.jac = function(par) ui,
            control.outer = list(trace = FALSE))[c("par", "value", "gradient")]
    } else {
        opt <- alabama::auglag(par = start, fn = ll, gr = gr,
            hin = function(par) ui %*% par - ci, 
            hin.jac = function(par) ui,
            control.outer = list(trace = FALSE))[c("par", "value", "gradient")]
    }
    if (length(eY$which) > 0) {
        gamma <- opt$par[-(1:ncol(eY$Y))]
        names(opt$par)[-(1:ncol(eY$Y))] <- paste0("gamma", 1:length(gamma))
    } else {
        gamma <- opt$par[-(1:ncol(iY$Yleft))]
        names(opt$par)[-(1:ncol(iY$Yleft))] <- paste0("gamma", 1:length(gamma))
    }
    Lambdat@x[] <- mapping(gamma)
    opt$G <- crossprod(Lambdat)[1:length(rt$cnms[[1]]),1:length(rt$cnms[[1]])]
    if (Hessian) opt$Hessian <- numDeriv::hessian(ll, opt$par)
    opt$loglik <- ll
    class(opt) <- "mtramc"
    opt
}

logLik.mtramc <- function(object, parm = NULL, ...) {
    if (!is.null(parm)) {
        ret <- -c(object$loglik(parm))
    } else {
        ret <- -c(object$value)
    }
    attr(ret, "df") <- length(coef(object))
    class(ret) <- "logLik"
    ret
}

coef.mtramc <- function(object, ...)
    object$par

.Marsaglia_1963 <- function(lower = rep(-Inf, nrow(sigma)), 
                            upper = rep(Inf, nrow(sigma)), 
                            mean = rep(0, nrow(sigma)), 
                            V = diag(2), 
                            grd = NULL,
                            do_qnorm = TRUE,
                            ...) {
 
    k <- nrow(V)
    l <- ncol(V)

    if (is.null(grd)) {
        stopifnot(do_qnorm)
        grd <- SparseGrid::createSparseGrid(type = "KPU", dimension = ncol(V), k = 10)
        ### Note: We expect t(nodes) below
        grd$nodes <- t(grd$nodes)
    }
    
    lower <- lower - mean
    upper <- upper - mean

    if (k == 1) {
        VVt <- base::tcrossprod(V)
        sd <- sqrt(base::diag(VVt) + 1)
        return(pnorm(upper / sd) - pnorm(lower / sd))
    }

    ### y = qnorm(x)
    inner <- function(y) {
        Vy <- V %*% y
        ### this needs ~ 75% of the total runtime
        #ret <- pnorm(upper - Vy) - pnorm(lower - Vy)
        ### ~ 3x speed-up
        ### ret <- .Call("pnormMRS", c(upper - Vy)) - .Call("pnormMRS", c(lower - Vy))
        #if (nrow(Vy) == 1) return(ret)
        #ret <- matrix(pmax(.Machine$double.eps, ret), nrow = nrow(Vy),
        #              ncol = ncol(Vy))
        #exp(colSums(log(ret)))
        .Call("R_inner", upper - Vy, lower - Vy)
    }

    if (do_qnorm) grd$nodes <- qnorm(grd$nodes)
    ev <- inner(grd$nodes)
    c(value = sum(grd$weights * ev))
}
