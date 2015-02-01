
as.bases <- function(object, ...)
    UseMethod("as.bases")

as.bases.formula <- function(object, remove_intercept = FALSE, 
                             data = NULL, 
                             ci = NULL, ui = NULL, ...) {

    varnames <- all.vars(object)

    ret <- function(data, ...) {
        data <- data[varnames]
        data <- as.data.frame(data)
        mf <- model.frame(object, data, ...)
        X <- model.matrix(object, data = mf, ...)
        if (remove_intercept)
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        if (is.null(ui)) ui <- Diagonal(ncol(X))
        if (is.null(ci)) ci <- rep(-Inf, ncol(X))
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        return(X)
    }
    if (!is.null(data)) {
        s <- lapply(data, function(v) {
            if (is.factor(v)) return(unique(v))
            range(v)
        })
    } else {
        s <- NULL
    }
    attr(ret, "varnames") <- varnames
    attr(ret, "support") <- s
    class(ret) <- c("bases", "basis", class(ret))
    ret
}
