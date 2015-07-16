
polynomial_basis <- function(coef, support = c(0, 1),
                             ui = NULL, ci = NULL,
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

    if (is.null(ui)) ui <- Diagonal(length(object))
    if (is.null(ci)) ci <- rep(-Inf, length(object))

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
        if (deriv < 0) X[] <- 0
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

orthopolynomial_basis <- function(order, support = c(0, 1),
                                  type = c("Legendre"),
                                  ui = c("none", "increasing", "decreasing"), # "cyclic", "zerointegral"),
                                  varname = NULL, ...) {

    type <- match.arg(type)
    object <- legendre.polynomials(order, ...)

    ui <- match.arg(ui)
    if (type == "Legendre") {
        B <- Bernstein_basis(order, ui = ui, varname = varname)
        constr <- get("constr", environment(B))
    } else {
        ui <- match.arg(ui)
        stopifnot(ui == "none")
        constr <- list(ui = Diagonal(order + 1), ui = rep(-Inf, order + 1))
    }

    basis <- function(data, deriv = 0L) {
        stopifnot(order > max(c(0, deriv)))
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

        ### map into [0, 1]
        x <- (x - support[1]) / diff(support)
        stopifnot(all(x >= 0 && x <= 1))

        ### Legendre is on [-1, 1] !
        X <- do.call("cbind", as.vector(lapply(dobject, predict, 2 * x - 1)))
        colnames(X) <- paste("L", 1:ncol(X), sep = "")
        if (deriv > 0)
            X <- X * (2 / diff(support)^deriv)
        if (deriv < 0) X[] <- 0
        if (type == "Legendre" && ui != "none") {
            L2B <- .Call("L2B", as.integer(order))
            constr$ui <- constr$ui %*% L2B
        }
        attr(X, "constraint") <- constr
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
    object(data, deriv = .deriv(variable.names(object), deriv))
