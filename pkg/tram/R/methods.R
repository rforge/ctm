
as.mlt.tram <- function(object) {
    cls <- which(class(object) == "mlt")
    class(object) <- class(object)[-(1:(cls - 1))]
    object
}    


model.matrix.tram <- function(object, baseline = FALSE, ...) 
{
    ret <- model.matrix(as.mlt(object), ...)
    if (baseline) return(ret)
    return(ret[, object$shiftcoef,,drop = FALSE])
}	

coef.tram <- function(object, baseline = FALSE, ...) 
{
    cf <- coef(as.mlt(object), ...)
    if (baseline) return(cf)
    return(cf[object$shiftcoef])
}
        
vcov.tram <- function(object, baseline = FALSE, ...) 
{
    if (is.null(object$cluster)) {
        ret <- vcov(as.mlt(object), ...)
    } else {
        ret <- sandwich::vcovCL(as.mlt(object), cluster = object$cluster)
    }
    if (baseline) return(ret)
    return(ret[object$shiftcoef, object$shiftcoef, drop = FALSE])
}

nobs.tram <- function(object, ...)
    NROW(object$data)

logLik.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
    logLik(as.mlt(object), parm = parm, ...)

Hessian.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
    Hessian(as.mlt(object), parm = parm, ...)

Gradient.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
   Gradient(as.mlt(object), parm = parm, ...)

estfun.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
    estfun(as.mlt(object), parm = parm, ...)

predict.tram <- function(object, newdata = object$data, 
    type = c("lp", "trafo", "distribution", "survivor", "density", 
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"), ...) {

    type <- match.arg(type)
    if (type == "lp")
        return(model.matrix(object, data = newdata) %*% coef(object, baseline = FALSE))
    predict(as.mlt(object), newdata = newdata, type = type, ...)
}

print.tram <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(coef(x, baseline = FALSE))
    ll <- logLik(x)
    cat("\nLog-Likelihood:\n ", ll, " (df = ", attr(ll, "df"), ")", sep = "")
    cat("\n\n")
    invisible(x)
}

summary.tram <- function(object, ...) {
    ret <- list(call = object$call,
                test = cftest(object, parm = names(coef(object, baseline = FALSE))),
                ll = logLik(object))
    class(ret) <- "summary.tram"
    ret
}

print.summary.tram <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    pq <- x$test$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    sig <- .Machine$double.eps
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
        P.values = TRUE, eps.Pvalue = sig)
    cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), ")", sep = "")
    cat("\n\n")
    invisible(x)
}
