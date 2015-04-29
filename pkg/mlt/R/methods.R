
coef.mlt <- function(object, fixed = TRUE, ...) {
    if (fixed) 
        return(object$coef)
    return(object$par)
}

"coef<-" <- function(object, value)
    UseMethod("coef<-")

"coef<-.mlt" <- function(object, value) {
    cf <- coef(object, fixed = TRUE)
    stopifnot(length(cf) == length(value))
    if (!is.null(names(value)))
        stopifnot(all.equal(names(cf), names(value)))
    object$par <- object$parm(value)
    object$coef <- value
    object
}

### use maxLik::hessian?
Hessian <- function(object, ...)
    UseMethod("Hessian")

Hessian.mlt <- function(object, ...)
    object$hessian(coef(object, fixed = FALSE))
    
vcov.mlt <- function(object, ...)
    solve(Hessian(object, ...))

logLik.mlt <- function(object, ...) {
    ret <- -object$loglik(coef(object, fixed = FALSE))
    attr(ret, "df") <- length(coef(object, fixed = FALSE))
    class(ret) <- "logLik"
    ret
}

estfun.mlt <- function(object, ...)
    -object$score(coef(object, fixed = FALSE))

mkgrid.mlt <- function(object, ...)
    mkgrid(object$model, ...)

mkgrid.model <- function(object, ...)
    mkgrid(object$model, ...)

variable.names.mlt <- function(object, ...)
    variable.names(object$model)
