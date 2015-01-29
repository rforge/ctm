
.mlt_fit <- function(response, data, bresponse, binter = NULL, bshift = NULL, distr) {
#, offsetfun, offset, trunc = c(-Inf, Inf)) {

    y <- data[[response]]
    tmpdata <- data
    if (!is.null(binter)) {
        by <- c(bresponse = bresponse, binter = binter)
    } else {
        by <- c(bresponse = bresponse)
    }
    if (!is.Surv(y)) {
        Y <- model.matrix(by, data = data)
        Yprime <- model.matrix(by, data = data, bresponse = list(deriv = 1))
        ui <- attr(Y, "constraint")$ui
        ci <- attr(Y, "constraint")$ci
        stopifnot(any(is.finite(ci)))
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        ci[ci == 0] <- sqrt(.Machine$double.eps)
        if (!is.null(bshift)) {
            X <- model.matrix(bshift, data = data) ### remove intercept
            if ("(Intercept)" %in% colnames(X))
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            Y <- cbind(Y, -X)
            Yprime <- cbind(Yprime, matrix(0, nrow = nrow(Yprime), ncol = ncol(X)))
            ui <- cbind(ui, matrix(0, nrow = nrow(ui), ncol = ncol(X)))
        }
        ll <- function(beta) -sum(.mlt_loglik(distr, Y, Yprime)(beta))
        sc <- function(beta) -colSums(.mlt_score(distr, Y, Yprime)(beta))
        z <- scale(y)
        theta <- coef(lm(z ~ Y - 1))
        ret <- constrOptim(theta = theta, f = ll, 
                           grad = sc, ui = ui,
                           ci = ci, hessian = TRUE)
    }
    ret$scores <- .mlt_score(distr, Y, Yprime)(ret$par)
    ret$by <- by
    ret$bshift <- bshift
    ret$response <- response
    ret$distr <- distr
    class(ret) <- "mlt"
    return(ret)
}

mlt <- .mlt_fit

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
        Y <- model.matrix(object$by, data = newdata)
        if (!is.null(object$bshift)) {
            X <- model.matrix(object$bshift, data = newdata) ### remove intercept
            if ("(Intercept)" %in% colnames(X))
                X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            Y <- cbind(Y, -X)
        }
        ret <- drop(Y %*% object$par)
        ret <- switch(type, 
                  "trafo" = ret,
                  "prob" = object$distr$p(ret))
        ret
    }

    ret
}
