
coef.mlt <- function(object, fixed = TRUE, ...) {
    if (fixed) 
        return(object$coef)
    return(object$par)
}

vcov.mlt <- function(object, ...)
    solve(object$optim(coef(object, fixed = FALSE), 
                       hessian = TRUE)$hessian)

logLik.mlt <- function(object, ...) {
    ret <- -object$loglik(coef(object, fixed = FALSE))
    attr(ret, "df") <- length(coef(object, fixed = FALSE))
    class(ret) <- "logLik"
    ret
}

estfun.mlt <- function(object, ...)
    -object$score(coef(object, fixed = FALSE))

predict.model <- function(object, newdata, coef, ...) {

    stopifnot(is.data.frame(newdata))
    lapply(1:nrow(newdata), function(i) {
        ret <- function(y, type = c("trafo", "prob", "density", 
                                    "quantile", "hazard", "cumhaz"), n = 1000, ...) { 

            type <- match.arg(type)

            tpdfun <- function(y, type = c("trafo", "prob", "density", "cumhaz")) {
                type <- match.arg(type)
                yt <- .type_of_response(y)
            
                data <- newdata[rep(i, length(y)),,drop = FALSE]
                data[[object$response]] <- y
                Y <- model.matrix(object, data = data, ...)
                trafo <- drop(Y %*% coef) ### <FIXME> offset??? </FIXME>
                ret <- switch(type, 
                    "trafo" = trafo,
                    "prob" = object$todistr$p(trafo),
                    "density" = {
                        if (yt %in% c("unordered", "ordered")) 
                            stop("not yet implemented")
                        deriv <- 1
                        names(deriv) <- object$response
                        Yprime <- model.matrix(object, data = data, deriv = deriv, ...)
                        object$todistr$d(trafo) * drop(Yprime %*% coef)
                    },
                    "cumhaz" = -log(1 - object$todistr$p(trafo))
                )
                return(ret)
            }
 
            if (type == "hazard")
                return(tpdfun(y, "density") / (1 - tpdfun(y, "prob")))
            if (type == "quantile") {
                probs <- y
                stopifnot(all(probs >= 0 & probs <= 1))
                y <- generate(object, n = n)[[object$response]]
                yt <- .type_of_response(y)
                p <- c(0, tpdfun(y, "prob") , 1)
                idx <- sapply(probs, function(prob) sum(p <= prob))
                if (yt %in% c("unordered", "ordered", "integer")) 
                    return(y[idx]) 
                ptmp <- cbind(p[idx], p[idx + 1])
                beta <- (ptmp[,2] - ptmp[,1]) / (y[idx + 1] - y[idx])
                alpha <- ptmp[,1] - beta * y[idx]
                return((probs - alpha) / beta)
            }
            return(tpdfun(y, type))
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

samplefrom.model <- function(object, newdata = NULL, coef, 
                             n = 1, ngrid = 10, ny = 100, interval = FALSE, ...) {

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
        if (interval)
            return(data.frame(left = y[idx], right = y[idx + 1]))
        ptmp <- cbind(p[cbind(1:nrow(p), idx)], p[cbind(1:nrow(p), idx + 1)])
        beta <- (ptmp[,2] - ptmp[,1]) / (y[idx + 1] - y[idx])
        alpha <- ptmp[,1] - beta * y[idx]
        return((u - alpha) / beta)
    }
    ret <- lapply(1:n, .sample)
    if (interval) {
        ret <- do.call("rbind", ret)
    } else { 
        ret <- unlist(ret)
    }
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
                     what = c("trafo", "prob"), ### density, quantile, ...
                     plotfun = plot, ...) {

    what <- match.arg(what)
    stopifnot(length(formula) == 2)
    formula[[3]] <- formula[[2]]
    formula[[2]] <- as.name("p")
    trafo <- predict(x$model$model, newdata = newdata, 
                     coef = coef(x, fixed = TRUE))
    p <- switch(what, 
        "trafo" = trafo,
        "prob" = x$model$todistr$p(trafo))
    nd <- expand.grid(newdata)
    nd$p <- p
    plotfun(formula, data = nd, ...)
}

generate.mlt <- function(object, ...)
    generate(object$model, ...)

generate.model <- function(object, ...)
    generate(object$model, ...)
