
trforest <- function(object, part, data, parm, weights, modelsplit = FALSE, 
             control = ctree_control(teststat = "quad",
                                     testtype = "Univ", mincriterion = 0, ...), 
             ...) {

    if (missing(parm)) parm <- 1:length(coef(object))

    trfo <- function(data, weights)
        estfun(update(object, weights = weights,
                      theta = coef(object)))[, parm, drop = FALSE]

    split <- function(x, response, weights, mb) {

        if (is.ordered(x)) x <- unclass(x)
        if (is.factor(x)) {
            if(nlevels(x) == 2) return(1:2)
            stop("not yet implemented")
        }
        ox <- order(x)
        w <- cumsum(weights[ox])
        ux <- x[ox][(w > mb) & (w < (sum(weights) - mb))]
        sux <- sort(unique(ux))
        suxl <- rev(sux[sux <= median(sux)])
        suxr <- sux[sux > median(sux)]

        lll <- numeric(length(suxl))
        lmod <- rmod <- object
        for (i in 1:length(suxl)) {
            lll[i] <- logLik(lmod <- update(object, weights = weights * (x <= suxl[i]),
                                        theta = coef(lmod))) +
                  logLik(rmod <- update(object, weights = weights * (x > suxl[i]),
                                        theta = coef(rmod)))
        }
        llr <- numeric(length(suxr))
        lmod <- rmod <- object
        for (i in 1:length(suxr)) {
            llr[i] <- logLik(lmod <- update(object, weights = weights * (x <= suxr[i]),
                                            theta = coef(lmod))) +
                      logLik(rmod <- update(object, weights = weights * (x > suxr[i]),
                                            theta = coef(rmod)))
        }
        xsplit <- c(suxl, suxr)[which.max(c(lll, llr))]
        xsplit
    }


    mf <- model.frame(part, data = data)
    stopifnot(all(!names(mf) %in% variable.names(object)))
    if (missing(weights)) weights <- rep(1, NROW(mf))
    tmp <- trfo(data, weights)

    mf$response <- 0
    fm <- response ~ x
    fm[[3L]] <- part[[2L]]

    if (modelsplit) control$splitfun <- split
    cf <- cforest(fm, data = mf, ytrafo = trfo, weights = weights, 
                  control = control, ...)

    cf$model <- object
    class(cf) <- c("trforest", class(cf))
    cf
}

predict.trforest <- function(object, newdata, K = 20, 
    type = c("node", "coef", "trafo", "distribution", 
             "survivor", "density", "logdensity", 
             "hazard", "loghazard", "cumhazard", "quantile"), FUN = NULL, ...) {

    class(object) <- class(object)[-1L]
    type <- match.arg(type)
    if (type == "node") return(predict(object, newdata = newdata, type = "node", ...))

    q <- mkgrid(object$model, n = K)[[object$model$response]]

    if (is.null(FUN)) {
        if (type == "coef") {
            FUN <- function(response, weights)
                coef(update(object$model, weights = weights))
        } else {
            FUN <- function(response, weights)
                  predict(update(object$model, weights = weights), q = q, type = type, ...)
        }
    }
    predict(object, newdata = newdata, FUN = FUN, ...)
}

logLik.trforest <- function(object, newdata, ...) {

    if (missing(newdata)) newdata <- object$model$data

    cf <- predict(object, newdata = newdata, type = "coef", ...)
    mod <- mlt(object$model$model, data = newdata, dofit = FALSE)

    ll <- rep(0, NROW(newdata))
    for (i in 1:NROW(newdata)) {
        w <- rep(0, NROW(newdata))
        w[i] <- 1
        ll[i] <- logLik(mod, parm = cf[i,], w = w)
    }
    ret <- sum(ll)
    attr(ret, "df") <- NA
    class(ret) <- "logLik"
    ret
}
