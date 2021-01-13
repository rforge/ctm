    
### simulate from model object with data newdata
simulate.ctm <- function(object, nsim = 1, seed = NULL, 
                         newdata, K = 50, q = NULL,
                         interpolate = TRUE, bysim = TRUE, ...) {

    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    y <- variable.names(object, "response")
    if (is.null(q))
        q <- mkgrid(object, n = K)[[y]]

    stopifnot(!missing(newdata))
    if (is.data.frame(newdata)) {
        p <- runif(nsim * NROW(newdata))
        ### basically compute quantiles for p; see qmlt
        prob <- predict(object, newdata = newdata, q = q, type = "distribution")
        ### unconditional: newdata is ignored in the presence of q
        if (!is.matrix(prob))
            prob <- matrix(prob, ncol = 1)[,rep(1, NROW(newdata)), drop = FALSE]
        prob <- t(prob[, rep(1:NCOL(prob), nsim),drop = FALSE])
        discrete <- !inherits(as.vars(object)[[y]],
                              "continuous_var")
        bounds <- bounds(object)[[y]]
        ret <- .p2q(prob, q, p, interpolate = interpolate, bounds = bounds,
                    discrete = discrete)
        if (nsim > 1) {
            if (bysim) {
                tmp <- vector(mode = "list", length = nsim)
                for (i in 1:nsim) {
                    idx <- 1:NROW(newdata) + (i - 1) * NROW(newdata)
                    if (is.data.frame(ret)) 
                        tmp[[i]] <- ret[idx,]
                    else 
                        tmp[[i]] <- ret[idx]
                }
            } else {
                tmp <- vector(mode = "list", length = NROW(newdata))
                for (i in 1:NROW(newdata)) {
                    idx <- NROW(newdata) * 0:(nsim - 1) + i
                    if (is.data.frame(ret)) 
                        tmp[[i]] <- ret[idx,]
                    else 
                        tmp[[i]] <- ret[idx]
                }
            }
            ret <- tmp
        }
    } else {
        stop("not yet implemented")
    }
    return(ret)
}

simulate.mlt <- function(object, nsim = 1, seed = NULL, newdata = object$data, bysim = TRUE, ...) {
    ctmobj <- object$model
    coef(ctmobj) <- coef(object)
    simulate(ctmobj, nsim = nsim, seed = seed, newdata = newdata, bysim = bysim, ...)
}

paraboot <- function(object, ...)
    UseMethod("paraboot")

paraboot.mlt <- function(object, parm = coef(object), B = 100, ...) {

    coefs <- matrix(0, nrow = B, ncol = length(parm))
    colnames(coefs) <- names(parm)
    logLR <- numeric(B)

    ndf <- object$data
    y <- variable.names(object, "response")

    for (b in 1:B) {
        ndf[[y]] <- simulate(object, ...)
        m <- mlt(model = object$model, data = ndf, theta = parm, weights = object$weights, 
                 scale = object$scale, offset = object$offset, check = FALSE, checkGrad = FALSE)
        coefs[b,] <- coef(m)
        logLR[b] <- -2 * (logLik(m, parm = parm) - logLik(m))
    }

    object$coefs <- coefs
    object$logLR <- logLR
    class(object) <- c("mlt_paraboot", class(object))
    object
}

sim <- function(object, nsim = 1, seed = NULL, 
                newdata, K = 50, q = NULL,
                bysim = TRUE, ...) {

    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    if (is.data.frame(newdata)) {
        U <- matrix(runif(nsim * NROW(newdata)), nrow = nsim, 
                    byrow = TRUE) ### be compatible with mlt < 1.2-1
        U <- pmax(U, sqrt(.Machine$double.eps))
        U <- pmin(U, 1 - sqrt(.Machine$double.eps))

    } else {
        stop("not yet implemented")
    }

    y <- variable.names(object, "response")
    if (is.null(q))
        q <- mkgrid(object, n = K)[[y]]

    ret <- .p2q2(object = object, newdata = newdata, p = U, q = q)

    if (nsim > 1) {
        if (bysim) {
            tmp <- vector(mode = "list", length = nsim)
            for (i in 1:nsim) {
                idx <- 1:NROW(newdata) + (i - 1) * NROW(newdata)
                if (is.data.frame(ret)) 
                    tmp[[i]] <- ret[idx,]
                else 
                    tmp[[i]] <- ret[idx]
            }
        } else {
            tmp <- vector(mode = "list", length = NROW(newdata))
            for (i in 1:NROW(newdata)) {
                idx <- NROW(newdata) * 0:(nsim - 1) + i
                if (is.data.frame(ret)) 
                    tmp[[i]] <- ret[idx,]
                else 
                    tmp[[i]] <- ret[idx]
            }
        }
        ret <- tmp
    }
    return(ret)
}

### numerically invert distribution function
.p2q2 <- function(object, newdata, p, q) {

    stopifnot(!missing(newdata))

    ### basically compute quantiles for p; see qmlt
    prob <- predict(object, newdata = newdata, q = q, type = "distribution")
    ### unconditional: newdata is ignored in the presence of q
    if (!is.matrix(prob))
        prob <- matrix(prob, ncol = 1)[,rep(1, NROW(newdata)), drop = FALSE]
    y <- variable.names(object, "response")
    bounds <- bounds(object)[[y]]

    discrete <- !inherits(as.vars(object)[[y]],
                          "continuous_var")
    if (discrete) {
        ### use "old" code
        prob <- t(prob[, rep(1:NCOL(prob), nrow(p)),drop = FALSE])
        prob <- cbind(0, prob, 1)
        i <- rowSums(prob < as.vector(t(p)))
        return(q[i])
    }

    ### use spline/approx to evaluate quantile function
    nleft <- nexact <- nright <- 
        matrix(NA, nrow = nrow(p), ncol = NROW(newdata))

    for (i in 1:NCOL(prob)) {
             
        pr <- prob[,i]
        pr0 <- which(pr < sqrt(.Machine$double.eps))
        pr1 <- which(pr > 1 - sqrt(.Machine$double.eps))
        rmi <- c()
        if (length(pr0) > 0L)
            rmi <- c(rmi, pr0[-length(pr0)])
        if (length(pr1) > 0L)
            rmi <- c(rmi, pr1[-1L])
        if (length(rmi) > 1L) {
            qq <- q[-rmi]
            pr <- pr[-rmi]
        } else {
            qq <- q
        }
        s <- spline(x = qq, y = pr, method = "hyman")
        ynew <- approx(x = s$y, y = s$x, xout = p[,i], 
                       yleft = -Inf, yright = Inf)$y
        nleft[ynew == -Inf, i] <- bounds[1L]
        nright[ynew == -Inf, i] <- min(q)
        nleft[ynew == Inf, i] <- max(q)
        nright[ynew == Inf, i] <- bounds[2L]
        ynew[!is.finite(ynew)] <- NA
        nexact[, i] <- ynew
    }

    if (all(is.na(nleft)) && all(is.na(nright)))
        return(as.vector(t(nexact)))
    return(R(as.vector(t(nexact)), 
             cleft = as.vector(t(nleft)), 
             cright = as.vector(t(nright))))
}
