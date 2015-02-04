
coef.mlt <- function(object, fixed = TRUE, ...) {
    if (fixed) 
        return(object$coef)
    return(object$par)
}

vcov.mlt <- function(object, ...)
    solve(object$optim(coef(object, fixed = FALSE), hessian = TRUE)$hessian)

logLik.mlt <- function(object, ...)
    -object$loglik(coef(object, fixed = FALSE))

estfun.mlt <- function(object, ...)
    -object$score(coef(object, fixed = FALSE))

predict.mlt <- function(object, newdata = NULL, 
                        ...) {
    ret <- function(y, type = c("trafo", "prob"),...) {
        type <- match.arg(type)
        if (is.null(newdata)) {
            newdata <- data.frame(y)
            colnames(newdata) <- object$response
        } else {
            idx <- 1:nrow(newdata)
            tmp <- expand.grid(y = y, idx = idx)
            newdata <- newdata[tmp$idx,,drop = FALSE]
            newdata[[object$response]] <- tmp$y
        }
        Y <- model.matrix(object$model, data = newdata)
        ret <- drop(Y %*% object$par)
        ret <- switch(type, 
                      "trafo" = ret,
                      "prob" = object$todistr$p(ret))
        ret
    }
    ret
}
