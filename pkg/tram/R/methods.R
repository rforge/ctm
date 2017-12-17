
as.mlt.tram <- function(object) {
    cls <- which(class(object) == "mlt")
    class(object) <- class(object)[-(1:(cls - 1))]
    object
}    

model.frame.tram <- function(formula, ...)
    formula$data

model.matrix.tram <- function(object, baseline = FALSE, ...) 
{
    ret <- model.matrix(as.mlt(object), ...)
    if (baseline) return(ret)
    if (is.null(object$shiftcoef)) return(NULL)
    return(ret[, object$shiftcoef,,drop = FALSE])
}	

coef.tram <- function(object, baseline = FALSE, ...) 
{
    cf <- coef(as.mlt(object), ...)
    if (baseline) return(cf)
    if (is.null(object$shiftcoef)) return(NULL)
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
    if (is.null(object$shiftcoef)) return(NULL)
    return(ret[object$shiftcoef, object$shiftcoef, drop = FALSE])
}

nobs.tram <- function(object, ...) {
    if (!is.null(object$weights)) 
        return(sum(object$weights != 0))
    return(NROW(object$data))
}

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
    cat("\n", x$tram, "\n")
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
    if (!is.null(object$LRtest)) {
        ret$LRstat <- object$LRtest["LRstat"]
        ret$df <- floor(object$LRtest["df"])
        ret$p.value <- pchisq(object$LRtest["LRstat"], 
                              df = object$LRtest["df"], lower.tail = FALSE)
    }
    class(ret) <- "summary.tram"
    ret
}

print.summary.tram <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\n", x$tram, "\n")
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
    if (!is.null(x$LRstat))
        cat("\nLikelihood-ratio Test: Chisq =", x$LRstat, "on",
            x$df, "degrees of freedom; p =", format.pval(x$p.value, digits = digits, ...))
    cat("\n\n")
    invisible(x)
}
