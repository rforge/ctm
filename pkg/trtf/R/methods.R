
.R2vec <- function(object) {
    if (!inherits(object, "response")) return(object)
    ex <- object$exact
    le <- object$cleft
    ri <- object$cright
    ex[is.na(ex)] <- 0
    le[is.na(le) | !is.finite(le)] <- 0
    ri[is.na(ri) | !is.finite(ri)] <- 0
    ex + (le + (ri - le) / 2)
}


coef.trafotree <- function(object, ...)
    object$coef

logLik.trafotree <- function(object, newdata, ...) {

    cf <- coef(object)
    if (missing(newdata)) {
        ret <- sum(object$logLik)
    } else {
        nd <- factor(predict(object, newdata = newdata, type = "node", ...))
        ### set up unfitted model with newdata
        mltargs <- object$mltargs
        mltargs$data <- newdata
        mltargs$dofit <- FALSE
        mltmod <- do.call("mlt", mltargs)
        weights <- rep(1, nrow(newdata))
        ll <- numeric(nlevels(nd))
        for (i in 1:nlevels(nd)) {
            w <- numeric(nrow(newdata))
            w[nd == levels(nd)[i]] <- 1
            ll[i] <- logLik(mltmod, parm = cf[levels(nd)[i],], w = w)
        }
        ret <- sum(ll)
    }
    attr(ret, "df") <- length(cf)
    class(ret) <- "logLik"
    ret
}

logLik.traforest <- function(object, newdata, OOB = FALSE, coef = NULL, ...) {

    if (is.null(coef)) {
        if (missing(newdata)) {
            cf <- predict(object, OOB = OOB, type = "coef")
        } else {
            cf <- predict(object, newdata = newdata, type = "coef")
        }
    } else {
        cf <- coef
    }
    if (missing(newdata)) {
        mltmod <- object$model
        newdata <- object$data
        weights <- object$fitted[["(weights)"]]
        if (is.null(weights)) weights <- rep(1, nrow(newdata))
    } else {
        ### set up unfitted model with newdata
        mltargs <- object$mltargs
        mltargs$data <- newdata
        mltargs$dofit <- FALSE
        mltmod <- do.call("mlt", mltargs)
        weights <- rep(1, nrow(newdata))
    }

    ret <- sum(sapply(1:nrow(newdata), function(i) {
        w <- numeric(length(weights))
        w[i] <- weights[i]
        logLik(mltmod, parm = cf[[i]], w = w)
    }))
    attr(ret, "df") <- NA
    class(ret) <- "logLik"
    ret
}

predict.trafotree <- function(object, newdata, K = 20, q = NULL,
    type = c("node", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    ...) {

    type <- match.arg(type)
    tmp <- object
    class(tmp) <- class(tmp)[-1L]
    if (missing(newdata)) {
        nf <- predict(tmp, type = "node")
    } else {
        nf <- predict(tmp, newdata = newdata, type = "node")
    }
    if (type == "node") return(nf)
    nf <- factor(nf)
    if (type == "coef") return(object$coef[nf,,drop = FALSE])

    mod <- object$model
    if (is.null(q))
        q <- mkgrid(mod, n = K)[[mod$response]]

    if (missing(newdata)) newdata <- data_party(object)

    ### <FIXME> need .R2vec??? </FIXME>
    pr <- .R2vec(predict(mod, newdata = newdata, q = q, type = type, ...))
    if (!is.matrix(pr))
        pr <- matrix(pr, nrow = NROW(pr), ncol = NROW(newdata))
    for (nd in levels(nf)) {
        i <- nf == nd
        coef(mod) <- object$coef[nd,]
        pr[,i] <- .R2vec(predict(mod, newdata = newdata[i,], q = q,
                                 type = type, ...))
    } 
    pr
}

predict.traforest <- function(object,  newdata, mnewdata = data.frame(1), K = 20, q = NULL,
    type = c("weights", "node", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    OOB = FALSE, simplify = FALSE, trace = FALSE, ...) {

    type <- match.arg(type)
    tmp <- object
    class(tmp) <- class(tmp)[-1L]

    ptype <- type
    if (!(ptype %in% c("weights", "node"))) ptype <- "weights"
    if (missing(newdata)) {
        ret <- predict(tmp, OOB = OOB, type = ptype,
                       simplify = TRUE)
    } else {
        ret <- predict(tmp, newdata = newdata, type = ptype,
                       simplify = TRUE)
    }
    if (type %in% c("weights", "node")) return(ret)

    mod <- object$model
    if (is.null(q))
        q <- mkgrid(mod, n = K)[[mod$response]]
    mltmod <- object$mltobj
    thetastart <- coef(mltmod$object)

    ans <- vector(mode = "list", length = ncol(ret))
    names(ans) <- colnames(ret)
    cf <- vector(mode = "list", length = ncol(ret))
    names(cf) <- colnames(ret)

    if (trace) pb <- txtProgressBar(style = 3)
    for (i in 1:ncol(ret)) {
        if (trace) setTxtProgressBar(pb, í / ncol(ret))
        w <- ret[,i]
        if (!is.null(mltmod$iy)) 
            w <- libcoin::ctabs(mltmod$iy, weights = w)[-1L]
        if (i > 1) {
            imin <- which.min(cs <- colSums((ret[, 1:(i - 1), drop = FALSE] - w)^2))
            thetastart <- cf[[imin]]
        }
        umod <- update(mltmod$object, theta = thetastart, weights = w)
        cf[[i]] <- coef(umod)
        if (type != "coef")
            ans[[i]] <- predict(umod, q = q, newdata = mnewdata, type = type)
    } 
    if (trace) close(pb)
    if (type == "coef") return(cf)
    return(ans)
}

simulate.traforest <- function(object, nsim = 1, seed = NULL, newdata, 
                               OOB = FALSE, coef = NULL, ...) {

    if (is.null(coef)) {
        if (missing(newdata)) {
            cf <- predict(object, type = "coef", OOB = OOB)
            newdata <- object$data
        } else {
            cf <- predict(object, newdata = newdata, type = "coef")
        }
    } else {
        cf <- coef
        newdata <- object$data
    }
    if (is.list(cf)) cf <- do.call("rbind", cf)
    if (nrow(cf) != nrow(newdata)) stop("coef and newdata don't match")

    mod <- object$model
    ret <- vector(mode = "list", length = nrow(cf))
    for (i in 1:nrow(cf)) {
        coef(mod) <- cf[i,]
        ret[[i]] <- simulate(mod, nsim = nsim, seed = seed, 
                             newdata = newdata[i,,drop = FALSE], bysim = FALSE, ...)[[1]]
    }
    ans <- vector(mode = "list", length = nsim)
    if (any(sapply(ret, function(x) inherits(x, "response")))) {
        ret <- lapply(ret, function(x) {
            if (inherits(x, "response")) return(x)
            R(x)
        })
        for (j in 1:nsim) {
            for (i in 1:nrow(cf))
                ans[[j]] <- rbind(ans[[j]], ret[[i]][j,])
        }
    } else {
        for (j in 1:nsim) {
            for (i in 1:nrow(cf))
                ans[[j]] <- rbind(ans[[j]], ret[[i]][j])
        }
    }
    ans
}

simulate.trafotree <- simulate.traforest
