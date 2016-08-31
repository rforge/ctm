
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

    if (missing(newdata)) {
        if (is.null(object$logLik))
            stop("use trafotree(..., refit = TRUE)")
        ret <- sum(object$logLik)
    } else {
        nd <- factor(predict(object, newdata = newdata, type = "node", ...))
        mod <- mlt(object$model, data = newdata, dofit = FALSE)
        ll <- rep(0, nlevels(nd))
        for (i in 1:nlevels(nd)) {
            w <- rep(0, NROW(newdata)) ### newdata is always unweighted
            w[nd == levels(nd)[i]] <- 1
            ll[i] <- logLik(mod, parm = coef(object)[levels(nd)[i],], w = w)
        }
        ret <- sum(ll)
    }
    attr(ret, "df") <- length(coef(object))
    class(ret) <- "logLik"
    ret
}

logLik.traforest <- function(object, newdata, OOB = FALSE, ...) {

    if (missing(newdata)) {
        cf <- predict(object, OOB = OOB, type = "coef")
    } else {
        cf <- predict(object, newdata = newdata, type = "coef")
    }
    mod <- object$model
    if (missing(newdata)) {
        newdata <- object$data
        weights <- rep(1, nrow(newdata))
    } else {
        weights <- object$fitted[["(weights)"]]
        if (is.null(weights)) weights <- rep(1, nrow(newdata))
    }
    mltmod <- mlt(mod, data = newdata, dofit = FALSE)

    ret <- sum(sapply(1:nrow(newdata), function(i) {
        w <- weights
        w[-i] <- 0
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

        mf <- object$data[, variable.names(mod), drop = FALSE]
        if (nmax < Inf) {
            bdr <- BDR::BDR(mf, total = TRUE, nmax = nmax, complete.cases.only = TRUE,
                            as.interval = mod$response)
            iy <- c(bdr)
            mf1 <- as.data.frame(bdr)
            attr(iy, "levels") <- 1:nrow(mf1)
            wi <- libcoin::ctabs(iy, weights = integer(0))[-1L]
        } else {
            mf1 <- mf   
            wi <- rep(1, NROW(mf))
        }

        mltmod <- mlt(mod, data = mf1, weights = wi, ...)
        thetastart <- coef(mltmod)

        FUN <- function(response, weights) {
            if (nmax < Inf)
                weights <- libcoin::ctabs(iy, weights = weights)[-1L]
            umod <- update(mltmod, theta = thetastart, weights = weights)
            if (type == "coef") return(coef(umod))
            return(predict(umod, q = q, newdata = mnewdata, type = type))
        }
        ptype <- "response"
    } 

    return(predict(tmp, newdata = newdata, OOB = OOB, FUN = FUN, type = ptype,
                   simplify = simplify))
}

simulate.traforest <- function(object, nsim = 1, seed, newdata, cf, ...) {

    if (missing(newdata)) {
        if (missing(cf))  
            cf <- predict(object, type = "coef", ...)
    } else {
        if (missing(cf))
            cf <- predict(object, newdata = newdata, type = "coef", ...)
    }
    if (is.list(cf)) cf <- do.call("rbind", cf)

    mod <- object$model
    # mod <- mlt(mod, data = newdata, dofit = FALSE)
    ret <- vector(mode = "list", length = nrow(cf))
    for (i in 1:nrow(cf)) {
        coef(mod) <- cf[i,]
        ret[[i]] <- simulate(mod, nsim = nsim, seed = seed, 
                             newdata = data.frame(1), ...)
    }
    ret
}

simulate.trafotree <- simulate.traforest
