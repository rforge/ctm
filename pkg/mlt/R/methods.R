
coef.mlt <- function(object, fixed = TRUE, ...) {
    if (fixed) 
        return(object$coef)
    return(object$par)
}

vcov.mlt <- function(object, ...)
    solve(object$optim(coef(object, fixed = FALSE), 
                       hessian = TRUE)$hessian)

logLik.mlt <- function(object, ...)
    -object$loglik(coef(object, fixed = FALSE))

estfun.mlt <- function(object, ...)
    -object$score(coef(object, fixed = FALSE))

predict.mlt <- function(object, newdata = NULL, ...) {

    if (is.null(newdata)) {
        if (length(varnames(object$model)) == 1)
            newdata <- data.frame(const = 1)
        else 
            newdata <- object$data
    }
    lapply(1:nrow(newdata), function(i) {
        ret <- function(y, type = c("trafo", "prob"), ...) {
            type <- match.arg(type)
            data <- newdata[rep(i, length(y)),,drop = FALSE]
            data[[object$response]] <- y
            Y <- model.matrix(object$model, data = data, ...)
            ret <- drop(Y %*% coef(object))
            switch(type, 
                   "trafo" = ret,
                   "prob" = object$todistr$p(ret))
        }
        environment(ret) <- new.env()
        assign("i", i, environment(ret)) ### i changes in parent.frame()!
        ret
    })
}

samplefrom <- function(object, ...)
    UseMethod("samplefrom")

samplefrom.mlt <- function(object, newdata = NULL, n = 1, ngrid = 100, 
                           ...) {

    response <- object$response
    if (is.null(newdata)) newdata <- object$data
    if (is.list(newdata)) {
        newdata[[response]] <- NULL
        newdata <- do.call("expand.grid", newdata)
    }
    ygrid <- generate(object$model, n = ngrid)[[response]]
    if (is.factor(ygrid)) ygrid <- ygrid[-length(ygrid)]
    pfun <- predict(object, newdata = newdata)

    idx <- do.call("c", lapply(pfun, function(f) {
        p <- f(ygrid, type = "prob")
        pOK <- p[!duplicated(p)]
        pOK <- pOK[pOK > 0 & pOK < 1]
        unclass(cut(runif(n), breaks = c(0, pOK, 1))) + sum(p < min(pOK))
    }))
    if (is.factor(ygrid)) {
        ret <- ygrid[idx]
    } else {
        ret <- data.frame(left = c(-Inf, ygrid, Inf)[idx],
                          right = c(-Inf, ygrid, Inf)[idx + 1])
    }
    newdata <- newdata[rep(1:nrow(newdata), rep(n, nrow(newdata))),,drop = FALSE]
    newdata[[object$response]] <- ret
    newdata
}

