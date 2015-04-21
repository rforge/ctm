
coef.mlt <- function(object, fixed = TRUE, ...) {
    if (fixed) 
        return(object$coef)
    return(object$par)
}

"coef<-" <- function(object, value, ...)
    UseMethod("coef<-")

"coef<-.mlt" <- function(object, value, ...) {
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
                y <- mkgrid(object, n = n)[[object$response]]
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
        if (length(variable.names(object)) == 1)
            newdata <- data.frame(const = 1)
        else 
            newdata <- object$data
    }
    predict(object$model, newdata = newdata, 
            coef = coef(object, fixed = TRUE), ...)
}

plot.mlt <- function(x, formula, newdata = mkgrid(x, n = 25), 
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

mkgrid.mlt <- function(object, ...)
    mkgrid(object$model, ...)

mkgrid.model <- function(object, ...)
    mkgrid(object$model, ...)
