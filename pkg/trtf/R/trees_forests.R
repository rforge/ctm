
trafotree <- function(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...) {

    mltargs$model <- object

    ctmfit <- function(formula, data, weights, cluster, ctrl, ...) {
        weights <- model.weights(data)
        if (is.null(weights)) weights <- integer(0)
        f <- Formula(formula)    
        mf <- model.frame(formula = f, data = data, na.action = na.pass)
        if (ctrl$nmax < Inf) {
            bdr <- inum::inum(mf, total = TRUE, nmax = ctrl$nmax, complete.cases.only = TRUE,
                              as.interval = object$response)
            iy <- c(bdr)
            mf1 <- as.data.frame(bdr)
            attr(iy, "levels") <- 1:nrow(mf1)
            w <- libcoin::ctabs(iy, weights = weights)[-1L]
        } else {
            mf1 <- mf
            w <- weights
            if (length(w) == 0) w <- rep(1, nrow(mf1))
            iy <- NULL
        }
        mltargs$data <- mf1  
        mltargs$weights <- w     
        object <- do.call("mlt", mltargs)
        thetastart <- coef(object)

        function(subset = NULL, info = NULL, newweights = NULL, model = FALSE, estfun = TRUE) {
            if (model) return(list(object = object, iy = iy))
            if (is.null(newweights)) {
                if (ctrl$nmax < Inf) {
                    w <- libcoin::ctabs(iy, weights = weights, subset = subset)[-1L]
                } else {
                    w[-subset] <- 0
                }
            } else {
                w <- newweights
            }
            if (!is.null(info$coef)) thetastart <- info$coef
            umod <- suppressWarnings(try(update(object, theta = thetastart, weights = w), silent = TRUE))
            if (inherits(umod, "try-error") || umod$convergence != 0) {
                umod <- suppressWarnings(try(update(object, weights = w), silent = TRUE))
                if (inherits(umod, "try-error") || umod$convergence != 0) {
                    mltargs$weights <- w
                    umod <- try(do.call("mlt", mltargs))
                }
            }
            if (inherits(umod, "try-error")) {
                return(list(estfun = matrix(0, nrow = nrow(mf1), ncol = length(w)),
                            iy = iy, converged = FALSE))
            }
            ret <- NULL
            if (estfun) { 
                ret <- estfun(umod)[, parm, drop = FALSE]
                ret <- ret / w ### we need UNWEIGHTED scores
                ret[w == 0,] <- 0
            }
            return(list(estfun = ret, index = iy, 
                        coef = coef(umod), logLik = logLik(umod), 
                        converged = isTRUE(all.equal(umod$convergence, 0))))
        }
    }

    ret <- ctree(..., ytrafo = ctmfit)
    ret$model <- object
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE)
    ret$mltargs <- mltargs

    nd <- predict(ret, type = "node")
    ret$models <- tapply(1:length(nd), factor(nd), function(i) 
        ret$trafo(i, estfun = FALSE)) ### note: trafo is (potentially) weighted
    ret$coef <- do.call("rbind", lapply(ret$models, function(x) x$coef))
    ret$logLik <- sapply(ret$models, function(x) x$logLik)

    class(ret) <- c("trafotree", class(ret))
    ret
}

traforest <- function(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...) {

    mltargs$model <- object
 
    ctmfit <- function(formula, data, weights, cluster, ctrl, ...) {
        weights <- model.weights(data)
        if (!is.null(weights)) 
            stop("trees in forests are unweighted")
        weights <- integer(0)
        f <- Formula(formula)    
        mf <- model.frame(formula = f, data = data, na.action = na.pass)
        if (ctrl$nmax < Inf) {
            bdr <- inum::inum(mf, total = TRUE, nmax = ctrl$nmax, complete.cases.only = TRUE,
                              as.interval = object$response)
            iy <- c(bdr)
            mf1 <- as.data.frame(bdr)
            attr(iy, "levels") <- 1:nrow(mf1)
            w <- libcoin::ctabs(iy, weights = weights)[-1L]
        } else {
            mf1 <- mf
            w <- weights
            if (length(w) == 0) w <- rep(1, nrow(mf1))
            iy <- NULL
        }
        mltargs$data <- mf1
        mltargs$weights <- w
        object <- do.call("mlt", mltargs)
        thetastart <- coef(object)

        function(subset, info = NULL, newweights = NULL, model = FALSE, estfun = TRUE) {
            if (model) return(list(object = object, iy = iy))
            if (is.null(newweights)) {
                if (ctrl$nmax < Inf) {
                    w <- libcoin::ctabs(iy, weights = weights, subset = subset)[-1L]
                } else {
                    w[-subset] <- 0
                }
            } else {
                w <- newweights
            }
            if (!is.null(info$coef)) thetastart <- info$coef
            umod <- suppressWarnings(try(update(object, theta = thetastart, weights = w), silent = TRUE))
            if (inherits(umod, "try-error") || umod$convergence != 0) {
                umod <- suppressWarnings(try(update(object, weights = w), silent = TRUE))
                if (inherits(umod, "try-error") || umod$convergence != 0) {
                    mltargs$weights <- w
                    umod <- try(do.call("mlt", mltargs))
                }
            }
            if (inherits(umod, "try-error")) {
                return(list(estfun = matrix(0, nrow = nrow(mf1), ncol = length(w)),
                            iy = iy, converged = FALSE))
            } 
            ret <- NULL
            if (estfun) {
                ret <- estfun(umod)[, parm, drop = FALSE]
                ret <- ret / w ### we need UNWEIGHTED scores
                ret[w == 0,] <- 0
            }
            return(list(estfun = ret, index = iy, 
                        coef = coef(umod), logLik = logLik(umod), 
                        converged = isTRUE(all.equal(umod$convergence, 0))))
        }
    }
    ret <- cforest(..., ytrafo = ctmfit)
    ret$model <- object
    ret$mltargs <- mltargs
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE)
    class(ret) <- c("traforest", class(ret))
    ret
}
