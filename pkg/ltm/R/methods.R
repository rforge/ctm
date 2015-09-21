
as.mlt.ltm <- function(object) {
    class(object) <- class(object)[-(1:2)]
    object
}

### is.ltm?

model.frame.ltm <- function(object, ...)
    object$data

model.matrix.ltm <- function(object, data = model.frame(object), lp_only = TRUE, ...) {
    if (!lp_only)
        return(model.matrix(as.mlt(object), data = data, ...))
    if (is.null(object$model$model$bshifting))
        return(NULL)
    model.matrix(object$model$model$bshifting, data = data, ...)
}

coef.ltm <- function(object, lp_only = TRUE, ...) {
    if (!lp_only)
        return(coef(as.mlt(object), ...))
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NULL)
    coef(as.mlt(object))[xn]
}

vcov.ltm <- function(object, lp_only = TRUE, ...) {
    if (!lp_only)
        return(vcov(as.mlt(object), ...))
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NULL)
    vcov(as.mlt(object))[xn, xn, drop = FALSE]
}

### do we really need this?
estfun.ltm <- function(object, lp_only = TRUE, ...) {
    if (!lp_only)
        return(estfun(as.mlt(object), ...))
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NULL)
    estfun(as.mlt(object))[, xn, drop = FALSE]
}

AIC.ltm <- function(object, ..., k = 2)
    AIC(as.mlt(object), ..., k = k)

logLik.ltm <- function(object, ...)
    logLik(as.mlt(object), ...)

paraboot.ltm <- function(object, ...)
    paraboot(as.mlt(object))

summary.ltm <- function(object, ...) {

    cf <- coef(as.mlt(object))
    cfltm <- coef(object)
    ret <- list(call = object$call, convergence = object$convergence)
    if (is.null(cfltm)) {
        ret$intercepts <- cf
        ret$coefficients <- NULL
    } else {
        ret$intercepts <- cf[!(names(cf) %in% names(cfltm))]
        x <- cftest(object, ...)
        pq <- x$test
        mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
        pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df == 
            0, "z", "t"), ")", sep = ""), greater = paste("Pr(>", 
            ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|", 
            ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
        colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df == 
            0, "z value", "t value"), pname)
        alt <- switch(x$alternative, two.sided = "==", less = ">=", 
            greater = "<=")
        rownames(mtests) <- paste(rownames(mtests), alt, x$rhs)
        ret$coefficients <- mtests
    }
    ### add chisq test
    ret$AIC <- AIC(object)
    ret$logLik <- logLik(object)
    class(ret) <- "summary.ltm"
    ret
}

print.summary.ltm <- function(x,  digits = max(3L, getOption("digits") - 3L), ...) {

    cat("\nCall:\n")
    print(x$call)
    if (x$convergence != 0L)
    cat("\nCould not estimate parameters; optimisation did not converge!\n")
    cat("\nCoefficients (response transformation):\n")
    print(x$intercepts)
    if (!is.null(x$coefficients)) {
        cat("\nCoefficients (linear predictor):\n")         
        printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE, 
                     P.values = TRUE, eps.Pvalue = .Machine$double.eps)
    }
    cat("\nAIC:", x$AIC)
    cat("\nLog-Likelihood: ", x$logLik, " (df = ", attr(x$logLik, "df"), ")", sep = "")
    cat("\n")
    invisible(x)
}

print.ltm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    if (x$convergence != 0L)
    cat("\nCould not estimate parameters; optimisation did not converge!\n")
    cat("\nCoefficients:\n")
    print(coef(x, lp_only = FALSE))
}
