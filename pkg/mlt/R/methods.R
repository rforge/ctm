
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
    object$coef[] <- value ### preserve names
    object
}

### use maxLik::hessian?
Hessian <- function(object, ...)
    UseMethod("Hessian")

Hessian.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    object$hessian(parm)
    
Gradient <- function(object, ...)
    UseMethod("Gradient")

Gradient.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    as.vector(colSums(estfun(object, parm = parm)))

vcov.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    solve(Hessian(object, parm = parm))

logLik.mlt <- function(object, parm = coef(object, fixed = FALSE), ...) {
    ret <- -object$loglik(parm)
    attr(ret, "df") <- length(coef(object, fixed = FALSE))
    class(ret) <- "logLik"
    ret
}

estfun.mlt <- function(object, parm = coef(object, fixed = FALSE), ...)
    -object$score(parm)

mkgrid.mlt <- function(object, n, ...)
    mkgrid(object$model, n = n, ...)

mkgrid.ctm <- function(object, n = n, bounds = NULL, ...) {
    ret <- mkgrid(object$model, n = n, ...)
    if (!is.null(bounds)) {
        stopifnot(length(bounds) == 1)
        stopifnot(names(bounds) %in% names(ret))
        m <- ret[[names(bounds)]]
        if (is.double(m)) { ### <FIXME> what about integers? </FIXME>
            b <- bounds[[1]]
            s <- range(m)
            if (is.finite(b[1])) {
                stopifnot(b[1] <= s[1])
                s[1] <- b[1]
            }
            if (is.finite(b[2])) {
                stopifnot(s[2] <= b[2])
                s[2] <- b[2]
            }
            ret[[names(bounds)]] <- seq(from = s[1], to = s[2], length.out = length(m))
        }
    }
    ret
}

variable.names.mlt <- function(object, ...)
    variable.names(object$model)

as.mlt <- function(object)
    UseMethod("as.mlt")
