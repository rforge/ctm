
as.basis.formula <- function(object, remove_intercept = FALSE, 
                             data = NULL, ci = NULL, ui = NULL, ...) {

    varnames <- all.vars(object)

    ret <- function(data) {
        data <- data[varnames]
        data <- as.data.frame(data)
        mf <- model.frame(object, data)
        X <- model.matrix(object, data = mf, ...)
        if (remove_intercept) {
            a <- attr(X, "assign")
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            attr(X, "assign") <- a[-1]
        }
        if (is.null(ui)) ui <- Diagonal(ncol(X))
        if (is.null(ci)) ci <- rep(-Inf, ncol(X))
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- c("(Intercept)", varnames)[attr(X, "assign") + 1]
        return(X)
    }
    if (!is.null(data)) {
        s <- lapply(data, function(v) {
            if (is.factor(v)) return(unique(v))
            if (length(v) > 0)
                return(range(v, na.rm = TRUE))
            return(NULL)
        })
    } else {
        s <- NULL
    }
    attr(ret, "varnames") <- varnames
    attr(ret, "support") <- s
    class(ret) <- c("basis", class(ret))
    ret
}
