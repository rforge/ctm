
### set-up Bernstein polynom basis functions
### http://en.wikipedia.org/wiki/Bernstein_polynomial
### http://dx.doi.org/10.1080/02664761003692423
Bernstein_basis <- function(order = 2, support = c(0, 1),
                            ui = c("none", "increasing", "decreasing", "cyclic", "zerointegral"), 
                            ci = 0, varname = NULL) {

    zeroint <- FALSE
    ui <- match.arg(ui, several.ok = TRUE)
    if (length(ui) > 1) 
        stopifnot(length(ui == 2) && "zerointegral" %in% ui)
    if ("zerointegral" %in% ui) {
        zeroint <- TRUE
        ui <- ifelse(length(ui) == 2, ui[ui != "zerointegral"], "none")
    }
    constr <- switch(ui,
        "none" = list(ui = Diagonal(order + 1), 
                      ci = rep(-Inf, order + 1)),
        "increasing" = list(ui = diff(Diagonal(order + 1), differences = 1), 
                            ci = rep(ci, order)),
        "decreasing" = list(ui = diff(Diagonal(order + 1), differences = 1) * -1,
                            ci = rep(ci, order)),
        "cyclic" = stop("not yet implemented"))

    basis <- function(data, deriv = 0L, integrate = FALSE) {
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        stopifnot(order > deriv)
        x <- (x - support[1]) / diff(support)
        stopifnot(all(x >= 0 && x <= 1))
        fun <- ifelse(isTRUE(integrate), pbeta, dbeta)
        X <- do.call("cbind", lapply(0:(order - deriv), 
            function(m) fun(x, shape1 = m + 1, 
                shape2 = (order - deriv) - m + 1) / ((order - deriv) + 1)))
        if (deriv > 0L) {
            fact <- prod(order:(order - deriv + 1)) * (1 / diff(support)^deriv)
            X <- X %*% diff(diag(order + 1), differences = deriv) * fact
        }
        colnames(X) <- 1:ncol(X)
        if (zeroint) {
            X <- X[, -ncol(X), drop = FALSE] - X[, ncol(X), drop = TRUE]
            constr$ui <- constr$ui[, -ncol(constr$ui),drop = FALSE]
        }
        attr(X, "constraint") <- constr
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        return(X)
    }

    s <- list(support)
    names(s) <- varname
    attr(basis, "support") <- s
    attr(basis, "varnames") <- varname

    class(basis) <- c("Bernstein_basis", "basis", class(basis))
    return(basis)
}

### evaluate model.matrix of Bernstein polynom
model.matrix.Bernstein_basis <- function(object, data,
                                         deriv = 0L, integrate = FALSE, ...)
    model.matrix.basis(object = object, data = data, 
                       deriv = deriv, integrate = integrate, ...)
