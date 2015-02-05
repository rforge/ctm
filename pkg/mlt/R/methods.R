
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

predict.model <- function(object, newdata, coef, ...) {

    stopifnot(is.data.frame(newdata))
    lapply(1:nrow(newdata), function(i) {
        ret <- function(y, type = c("trafo", "prob"), ...) { ### density, quantile, hazard, cumhaz
            type <- match.arg(type)
            data <- newdata[rep(i, length(y)),,drop = FALSE]
            data[[object$response]] <- y
            Y <- model.matrix(object, data = data, ...)
            ret <- drop(Y %*% coef)
            switch(type, 
                   "trafo" = ret,
                   "prob" = object$todistr$p(ret))
        }
        environment(ret) <- new.env()
        assign("i", i, environment(ret)) ### i changes in parent.frame()!
        ret
    })
}

predict.mlt <- function(object, newdata = NULL, ...) {

    if (is.null(newdata)) {
        if (length(varnames(object)) == 1)
            newdata <- data.frame(const = 1)
        else 
            newdata <- object$data
    }
    predict(object$model, newdata = newdata, 
            coef = coef(object, fixed = TRUE), ...)
}


samplefrom <- function(object, ...)
    UseMethod("samplefrom")

samplefrom.model <- function(object, newdata = NULL, coef, n = 1, ngrid = 10, ny = 100, ...) {

    response <- object$response

    if (is.null(newdata)) {
        simdata <- generate(object$model, n = ngrid)
        if (storage.mode(y) == "double")
            y <- seq(from = min(y), to = max(y), length = ny)
        if (is.factor(y)) y <- sort(y)[-nlevels(y)]
        simdata[[response]] <- y
        newdata <- simdata
        newdata[[response]] <- NULL
        newdata <- expand.grid(newdata)
    } else {
        y <- generate(object$model, n = ny)[[response]]
        if (is.factor(y)) y <- sort(y)[-nlevels(y)]
        simdata <- expand.grid(y = y, id = 1:nrow(newdata))
        simdata <- cbind(simdata, newdata[simdata$id,,drop = FALSE])
        simdata <- as.data.frame(simdata)
    }

    p <- object$todistr$p(predict(object$model, newdata = simdata, coef = coef, ...))
    ny <- length(y)
    nobs <- nrow(newdata)
    p <- matrix(p, nrow = nobs, ncol = ny, byrow = TRUE)
    p <- cbind(0, p, 1)

    .sample <- function(i) {
        idx <- rowSums(p <= (u <- runif(nobs)))
        if (storage.mode(y) != "double")
            return(y[idx])
        ptmp <- cbind(p[cbind(1:nrow(p), idx)], p[cbind(1:nrow(p), idx + 1)])
        beta <- (ptmp[,2] - ptmp[,1]) / (y[idx + 1] - y[idx])
        alpha <- ptmp[,1] - beta * y[idx]
        return((u - alpha) / beta)
    }
    ret <- unlist(lapply(1:n, .sample))
    newdata <- newdata[rep(1:nrow(newdata), n),,drop = FALSE]
    newdata[[response]] <- ret
    return(newdata)
}

samplefrom.mlt <- function(object, newdata = NULL, n = 1, ...) {

    if (is.null(newdata)) newdata <- object$data
    samplefrom(object$model, newdata = newdata, 
               coef = coef(object, fixed = TRUE), n = n, ...)
}

plot.mlt <- function(x, formula, newdata = generate(x, n = 25), 
                     type = c("trafo", "prob"), plotfun = plot, ...) {

    type <- match.arg(type)

    p <- predict(x$model$model, newdata = newdata, 
                 coef = coef(x, fixed = TRUE))
    if (type == "prob") p <- x$model$todist$p(p)
    nd <- expand.grid(newdata)
    nd$p <- p
    plotfun(formula, data = nd, ...)
}

generate.mlt <- function(object, ...)
    generate(object$model, ...)

generate.model <- function(object, ...)
    generate(object$model, ...)
