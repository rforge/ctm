
model.frame.ltm <- function(object, ...)
    object$data

model.matrix.ltm <- function(object, data = model.frame(object), all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(model.matrix(object, data = data, ...))
    }
    if (is.null(object$model$model$bshifting))
        return(NA)
    model.matrix(object$model$model$bshifting, data = data, ...)
}

coef.ltm <- function(object, all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(coef(object, ...))
    }      
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    class(object) <- class(object)[-1]
    coef(object)[xn]
}

vcov.ltm <- function(object, all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(vcov(object, ...))
    }      
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    class(object) <- class(object)[-1]
    vcov(object)[xn, xn, drop = FALSE]
}

estfun.ltm <- function(object, all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(estfun(object, ...))
    }      
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    class(object) <- class(object)[-1]
    estfun(object)[, xn, drop = FALSE]
}

AIC.ltm <- function(object, ..., k = 2) {
    class(object) <- class(object)[-1]
    AIC(object, ..., k = k)
}

logLik.ltm <- function(object, ...) {
    class(object) <- class(object)[-1]
    logLik(object, ...)
}

print.ltm <- function(x, ...) 
    print(cftest(x))

paraboot.ltm <- function(object, ...) {
    class(object) <- class(object)[-1]
    paraboot(object)
}
