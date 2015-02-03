
coef.mlt <- function(object, ...)
    object$par

vcov.mlt <- function(object, ...)
    solve(object$optim(coef(object), hessian = TRUE)$hessian)

logLik.mlt <- function(object, ...)
    object$loglik(coef(object))

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
