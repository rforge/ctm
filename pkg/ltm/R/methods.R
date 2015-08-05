
as.mlt.ltm <- function(object) {
    class(object) <- class(object)[-(1:2)]
    object
}

model.frame.ltm <- function(object, ...)
    object$data

model.matrix.ltm <- function(object, data = model.frame(object), all = FALSE, ...) {
    if (all)
        return(model.matrix(as.mlt(object), data = data, ...))
    if (is.null(object$model$model$bshifting))
        return(NA)
    model.matrix(object$model$model$bshifting, data = data, ...)
}

coef.ltm <- function(object, all = FALSE, ...) {
    if (all)
        return(coef(as.mlt(object), ...))
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    coef(as.mlt(object))[xn]
}

vcov.ltm <- function(object, all = FALSE, ...) {
    if (all)
        return(vcov(as.mlt(object), ...))
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    vcov(as.mlt(object))[xn, xn, drop = FALSE]
}

### do we really need this?
estfun.ltm <- function(object, all = FALSE, ...) {
    if (all)
        return(estfun(as.mlt(object), ...))
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    estfun(as.mlt(object))[, xn, drop = FALSE]
}

AIC.ltm <- function(object, ..., k = 2)
    AIC(as.mlt(object), ..., k = k)

logLik.ltm <- function(object, ...)
    logLik(as.mlt(object), ...)

print.ltm <- function(x, ...) {
    if (!is.na(coef(x)))
        print(cftest(x))
}

paraboot.ltm <- function(object, ...)
    paraboot(as.mlt(object))
