
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
            stopifnot(!("(Intercept)" %in% colnames(X)))
            Y <- cbind(Y, X)
            Yprime <- cbind(Yprime, X)
            ui <- cbind(ui, matrix(0, nrow = nrow(ui), ncol = ncol(X)))
        }
        ll <- function(beta) -sum(.mlt_loglik(distr, Y, Yprime))
        sc <- function(beta) -colSums(.mlt_score(distr, Y, Yprime))
        z <- scale(y)
        theta <- coef(lm(z ~ Y - 1))
        ret <- constrOptim(theta = theta, f = ll, 
                           grad = sc, ui = ui,
                           ci = ci, hessian = TRUE)
    }
    ret$scores <- .mlt_score(distr, Y, Yprime)(beta)
    return(ret)
}

mlt <- .mlt_fit

