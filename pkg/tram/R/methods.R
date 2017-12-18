
as.mlt.tram <- function(object) {
    cls <- which(class(object) == "mlt")
    class(object) <- class(object)[-(1:(cls - 1))]
    object
}    

model.frame.tram <- function(formula, ...)
    formula$data

model.matrix.tram <- function(object, with_baseline = FALSE, ...) 
{
    ret <- model.matrix(as.mlt(object), ...)
    if (with_baseline) return(ret)
    if (is.null(object$shiftcoef)) return(NULL)
    return(ret[, object$shiftcoef,,drop = FALSE])
}	

coef.tram <- function(object, with_baseline = FALSE, ...) 
{
    cf <- coef(as.mlt(object), ...)
    if (with_baseline) return(cf)
    if (is.null(object$shiftcoef)) return(NULL)
    return(cf[object$shiftcoef])
}

coef.Lm <- function(object, as.lm = FALSE, ...) {

    class(object) <- class(object)[-1L]
    if (!as.lm)
        return(coef(object, ...))

    if (!is.null(object$stratacoef))
        stop("cannot compute scaled coefficients with strata")

    cf <- coef(object, with_baseline = TRUE, ...)
    cfx <- coef(object, with_baseline = FALSE, ...)
    cfs <- cf[!(names(cf) %in% names(cfx))]
    sd <- 1 / cfs[names(cfs) != "(Intercept)"]

    if (is.null(object$shiftcoef))
        return(-cfs["(Intercept)"] * sd)
    return(c(-cfs["(Intercept)"], cfx) * sd)
}

coef.Survreg <- function(object, as.survreg = FALSE, ...)
    coef.Lm(object, as.lm = as.survreg, ...)
        
vcov.tram <- function(object, with_baseline = FALSE, ...) 
{
    if (is.null(object$cluster)) {
        ret <- vcov(as.mlt(object), ...)
    } else {
        ret <- sandwich::vcovCL(as.mlt(object), cluster = object$cluster)
    }
    if (with_baseline) return(ret)
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
        return(model.matrix(object, data = newdata) %*% coef(object, with_baseline = FALSE))
    predict(as.mlt(object), newdata = newdata, type = type, ...)
}

print.tram <- function(x, ...) {
    cat("\n", x$tram, "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(coef(x, with_baseline = FALSE))
    ll <- logLik(x)
    cat("\nLog-Likelihood:\n ", ll, " (df = ", attr(ll, "df"), ")", sep = "")
    cat("\n\n")
    invisible(x)
}

summary.tram <- function(object, ...) {
    ret <- list(call = object$call,
                test = cftest(object, parm = names(coef(object, with_baseline = FALSE))),
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
