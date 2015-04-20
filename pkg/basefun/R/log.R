
log_basis <- function(support = c(.Machine$double.eps, Inf),
                      ui = Diagonal(1), ci = 0, varname = NULL) {

    basis <- function(data, deriv = 0L) {
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        if (deriv == 0) {
            X <- matrix(log(pmax(x, .Machine$double.eps)), ncol = 1)
        } else if (deriv == 1) {
            X <- matrix(1 / x, ncol = 1)
        } else if (deriv > 1) {
            stop("not yet implemented")
        } else if (deriv < 0) {
            X <- matrix(0, ncol = 1, nrow = length(x))
        }
        colnames(X) <- paste("log(", varname, ")", sep = "")
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- varname
        X
    }

    s <- list(support)
    names(s) <- varname
    attr(basis, "support") <- s
    attr(basis, "varnames") <- varname
    attr(basis, "intercept") <- FALSE

    class(basis) <- c("log_basis", "basis", class(basis))
    return(basis)
}

model.matrix.log_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = .deriv(varnames(object), deriv))
