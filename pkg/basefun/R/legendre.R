
Legendre_basis <- function(order, support = c(0, 1),
                           ui = c("none", "increasing", "decreasing", "cyclic"), # "zerointegral"),
                           varname = NULL, ...) {

    object <- legendre.polynomials(order, ...)

    ui <- match.arg(ui)
    B <- Bernstein_basis(order, ui = ui, varname = varname)
    constr <- get("constr", environment(B))

    basis <- function(data, deriv = 0L, integrate = FALSE) {

        if (integrate) stop("Integration not yet implemented")

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
                dobject <- sapply(dobject, deriv)
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
        if (ui != "none")
            constr$ui <- constr$ui %*% L2B(order)
        attr(X, "constraint") <- constr
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        X
    }

    s <- list(support)
    names(s) <- varname
    attr(basis, "support") <- s
    attr(basis, "varnames") <- varname
    attr(basis, "intercept") <- TRUE

    class(basis) <- c("Legendre_basis", "Bernstein_basis", "basis", class(basis))
    return(basis)
}
