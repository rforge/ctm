
trforest <- function(object, part, data, parm, weights, modelsplit = FALSE, 
             control = ctree_control(teststat = "quad",
                                     testtype = "Univ", mincriterion = 0), 
             ...) {

    ret <- .trparty(object = object, part = part, data = data,
                    weights = weights, modelsplit = modelsplit,
                    control = control, FUN = cforest, ...)
    class(ret) <- c("trforest", class(ret))
    ret
}

predict.trforest <- function(object, newdata, K = 20, 
    type = c("node", "coef", "trafo", "distribution", 
             "survivor", "density", "logdensity", 
             "hazard", "loghazard", "cumhazard", "quantile"), 
    OOB = FALSE, FUN = NULL, ...) {

    class(object) <- class(object)[-1L]
    type <- match.arg(type)
    if (type == "node") 
        return(predict(object, newdata = newdata, type = "node", 
                       OOB = OOB))

    q <- mkgrid(object$model, n = K)[[object$model$response]]

    if (is.null(FUN)) {
        if (type == "coef") {
            FUN <- function(response, weights)
                coef(update(object$model, weights = weights))
        } else {
            FUN <- function(response, weights)
                predict(update(object$model, weights = weights), q = q, 
                        type = type, ...)
        }
    }
    predict(object, newdata = newdata, FUN = FUN, OOB = TRUE)
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
