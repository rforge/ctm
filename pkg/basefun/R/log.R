
log_basis <- function(support = c(.Machine$double.eps, Inf),
                      ui = 1,
                      ci = 0,
                      varname = NULL, ...) {

    basis <- function(data, deriv = 0L) {
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- 1
            x <- data[[varname]]
        }
        if (deriv == 0) {
            ret <- matrix(log(x), ncol = 1)
        } else if (deriv == 1) {
            ret <- matrix(1 / x, ncol = 1)
        } else {
            ret <- matrix(NA, ncol = 1)
        }
        attr(ret, "constraint") <- list(ui = ui, ci = ci)
        ret
    }

    s <- list(support)
    names(s) <- varname
    attr(basis, "support") <- s
    attr(basis, "varnames") <- varname

    class(basis) <- c("log_basis", "basis", class(basis))
    return(basis)
}
