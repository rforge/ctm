
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
    OOB = FALSE, simplify = FALSE, nmax = Inf, ...) {

    type <- match.arg(type)
    tmp <- object
    class(tmp) <- class(tmp)[-1L]

    FUN <- NULL
    ptype <- type
    if (!(type %in% c("weights", "node"))) {
        mod <- object$model
        if (is.null(q))
            q <- mkgrid(mod, n = K)[[mod$response]]

        mltmod <- object$mltobj
        thetastart <- coef(mltmod$object)

        FUN <- function(response, weights) {
            if (!is.null(mltmod$iy)) 
                weights <- libcoin::ctabs(mltmod$iy, weights = weights)[-1L]
            umod <- update(mltmod$object, theta = thetastart, weights = weights)
            if (type == "coef") return(coef(umod))
            return(predict(umod, q = q, newdata = mnewdata, type = type))
        }
        ptype <- "response"
    } 

    if (missing(newdata))
        return(predict(tmp, OOB = OOB, FUN = FUN, type = ptype,
                       simplify = simplify))
    return(predict(tmp, newdata = newdata, FUN = FUN, type = ptype,
                   simplify = simplify))
}

simulate.traforest <- function(object, nsim = 1, seed = NULL, newdata, 
                               mnewdata = data.frame(1), coef = NULL, ...) {

    if (is.null(coef)) {
        if (missing(newdata)) {
            cf <- predict(object, type = "coef", ...)
        } else {
            cf <- predict(object, newdata = newdata, type = "coef", ...)
        }
    } else {
        cf <- coef
    }
    if (is.list(cf)) cf <- do.call("rbind", cf)

    mod <- object$model
    ret <- vector(mode = "list", length = nrow(cf))
    for (i in 1:nrow(cf)) {
        coef(mod) <- cf[i,]
        ret[[i]] <- simulate(mod, nsim = nsim, seed = seed, 
                             newdata = mnewdata, ...)
    }
    ret
}

simulate.trafotree <- simulate.traforest
