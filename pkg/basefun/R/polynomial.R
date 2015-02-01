

polynomial_basis <- function(order, support = c(0, 1),
                             ui = Diagonal(length(object)), 
                             ci = rep(-Inf, length(object)),
                             varname = NULL, ...) {

    object <- polynomial(rep(1, order))

    basis <- function(data, deriv = 0L) {
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- 1
            x <- data[[varname]]
        }
        tmp <- object
        if (deriv > 0) {
            for (i in 1:deriv)
                tmp <- deriv(tmp)
            zero <- matrix(0, nrow = length(x), ncol = deriv)
        }
        ret <- do.call("cbind", lapply(1:length(tmp), function(i) x^(i - 1) * coef(tmp)[i]))
        if (deriv > 0)
            ret <- cbind(zero, ret)
        attr(ret, "constraint") <- list(ui = ui, ci = ci)
        ret
    }

    s <- list(support)
    names(s) <- varname
    attr(basis, "support") <- s
    attr(basis, "varnames") <- varname

    class(basis) <- c("polynomial_basis", "basis", class(basis))
    return(basis)
}
