
trforest <- function(object, part, data, parm, weights, modelsplit = FALSE, 
             control = ctree_control(teststat = "quad",
                                     testtype = "Univ", mincriterion = 0), 
             ...) {

    ret <- .trparty(object = object, part = part, data = data, parm = parm,
                    weights = weights, modelsplit = modelsplit,
                    control = control, FUN = cforest, ...)
    class(ret) <- c("trforest", class(ret))
    ret
}

predict.trforest <- function(object, newdata, K = 20, 
    type = c("node", "weights", "coef", "trafo", "distribution", 
             "survivor", "density", "logdensity", 
             "hazard", "loghazard", "cumhazard", "quantile"), 
    OOB = FALSE, FUN = NULL, cluster = TRUE, ...) {

    class(object) <- class(object)[-1L]
    type <- match.arg(type)
    if (missing(newdata)) newdata <- NULL ### FIXME
    if (type %in% c("node", "weights")) 
        return(predict(object, newdata = newdata, type = type, 
                       OOB = OOB))

    w <- predict(object, newdata = newdata, type = "weights", OOB = OOB)
    q <- mkgrid(object$model, n = K)[[object$model$response]]

    if (is.null(FUN)) {
        if (type == "coef") {
            FUN <- function(theta) function(response, weights)
                coef(update(object$model, weights = weights, theta = theta))
        } else {
            FUN <- function(theta) function(response, weights)
                predict(update(object$model, weights = weights, theta = theta), q = q, 
                        type = type, ...)
        }
    }

    if (cluster) {
        cl <- factor(kmeans(t(w), center = floor(sqrt(NCOL(w))))$cluster)
    } else {
        cl <- factor(rep(1, NCOL(w)))
    }
    ret <- FUN(coef(object$model))(NA, weights(object$model))
    ret <- matrix(NA, nrow = length(ret), ncol = NCOL(w))
    for (n in levels(cl)) {
        i <- (cl == n)
        theta <- coef(update(object$model, weights = rowMeans(w[,i,drop = FALSE])))
        ret[,i] <- apply(w[,i,drop = FALSE], 2, function(w) FUN(theta)(NA, w))
    }
    ret
}

logLik.trforest <- function(object, newdata, cf, ...) {

    if (missing(newdata)) {
        if (missing(cf)) 
            cf <- predict(object, type = "coef", ...)
        mod <- object$model
    } else {
        if (missing(cf)) 
            cf <- predict(object, newdata = newdata, type = "coef", ...)
        mod <- mlt(object$model$model, data = newdata, dofit = FALSE)
    }

    ll <- rep(0, ncol(cf))
    for (i in 1:ncol(cf)) {
        w <- rep(0, ncol(cf))
        w[i] <- 1
        ll[i] <- logLik(mod, parm = cf[,i], w = w)
    }
    ret <- sum(ll)
    attr(ret, "df") <- NA
    class(ret) <- "logLik"
    ret
}
