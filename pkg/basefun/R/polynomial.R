
polynomial_basis <- function(coef, support = c(0, 1),
                             ui = Diagonal(length(object)), 
                             ci = rep(-Inf, length(object)),
                             varname = NULL, orthogonal = NULL) {

    if (!is.null(orthogonal)) {
        stopifnot(is.atomic(orthogonal))
        object <- poly.orth(orthogonal, degree = length(coef))
    } else {
        stopifnot(all(coef %in% c(0, 1)))
        object <- do.call("polylist", lapply(1:length(coef), function(i) {
            cf <- coef[1:i]
            cf[1:i < i] <- 0
            return(polynomial(cf))
        }))
    }

    basis <- function(data, deriv = 0L) {
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        dobject <- object
        if (deriv > 0) {
            for (i in 1:deriv)
                dobject <- deriv(dobject)
        }
        X <- sapply(dobject, predict, x)
        if (!is.matrix(X)) X <- matrix(X, nrow = 1)
        cn <- c("(Intercept)", varname)
        if (ncol(X) > 2)
            cn <- c(cn, paste(varname, "^", 2:(ncol(X) - 1), sep = ""))
        colnames(X) <- cn[1:ncol(X)]
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        X
    }

    s <- list(support)
    names(s) <- varname
    attr(basis, "support") <- s
    attr(basis, "varnames") <- varname
    attr(basis, "intercept") <- TRUE

    class(basis) <- c("polynomial_basis", "basis", class(basis))
    return(basis)
}

model.matrix.polynomial_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = deriv, ...)
