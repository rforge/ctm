
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
        ret <- sum(object$logLik)
    } else {
        nd <- factor(predict(object, newdata = newdata, type = "node", ...))
        mod <- mlt(object$model, data = newdata, dofit = FALSE)
        ll <- rep(0, nlevels(nd))
        for (i in 1:nlevels(nd)) {
            w <- rep(0, NROW(newdata))
            w[nd == levels(nd)[i]] <- 1
            ll[i] <- logLik(mod, parm = coef(object)[levels(nd)[i],], w = w)
        }
        ret <- sum(ll)
    }
    attr(ret, "df") <- length(coef(object))
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

predict.traforest <- function(object,  newdata, K = 20, q = NULL,
    type = c("weights", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    OOB = FALSE, nmax = Inf,
    ...) {

    type <- match.arg(type)
    tmp <- object
    class(tmp) <- class(tmp)[-1L]

    if (missing(newdata)) {
        w <- predict(tmp, OOB = OOB, type = "weights")
    } else {
        w <- predict(tmp, newdata = newdata, type = "weights")
    }

    if (type == "weights") return(w)

    mod <- object$model
    if (is.null(q))
        q <- mkgrid(mod, n = K)[[mod$response]]

    if (missing(newdata)) newdata <- object$data

    mf <- object$data
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

    mltmod <- mlt(mod, data = mf1, weights = wi, maxit = 10000)
    thetastart <- coef(mltmod)

    lapply(1:ncol(w), function(i) {
        if (nmax < Inf) {
            wi <- libcoin::ctabs(iy, weights = w[,i])[-1L]
        } else {
            wi <- w[,i]
        }      
        umod <- update(mltmod, theta = thetastart, weights = wi)
        if (type == "coef") return(coef(umod))
        return(predict(umod, q = q, newdata = newdata, type = type))
    })
}

simulate.traforest <- function(object, nsim = 1, newdata, cf, ...) {

    if (missing(newdata)) {
        if (missing(cf))  
            cf <- predict(object, type = "coef", ...)
    } else {
        if (missing(cf))
            cf <- predict(object, newdata = newdata, type = "coef", ...)
    }
    if (is.list(cf)) cf <- do.call("rbind", cf)

    mod <- object$model
    # mod <- mlt(mod, data = newdata, doFit = FALSE)
    ret <- vector(mode = "list", length = nrow(cf))
    for (i in 1:nrow(cf)) {
        coef(mod) <- cf[i,]
        ret[[i]] <- simulate(mod, nsim = nsim, newdata = data.frame(1), ...)
    }
    ret
}

simulate.trafotree <- simulate.traforest
