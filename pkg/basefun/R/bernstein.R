
.B <- function(j, n, integrate = FALSE) {
    function(x) {
        if (j < 0 || (n - j) < 0)
            return(rep(0, length(x)))
        if (integrate)
            return(pbeta(x, shape1 = j + 1, shape2 = n - j + 1) / (n + 1))
        return(dbeta(x, shape1 = j + 1, shape2 = n - j + 1) / (n + 1))
    }
}

.Bx <- function(x, j, n, deriv = 0L, integrate = FALSE) {
    if (deriv == 0L)
        return(.B(j, n, integrate = integrate)(x))
    return(n * (.Bx(x, j - 1, n - 1, deriv = deriv - 1L, integrate = integrate) -
                .Bx(x, j,     n - 1, deriv = deriv - 1L, integrate = integrate)))
}

.oBx <- function(x, j, n, deriv = 0L, integrate = FALSE) {
    ret <- 0
    for (k in 0:j)
        ret <- ret + (-1)^k * choose(2 * n + 1 - k, j - k) *
               choose(j, k) * .Bx(x, j - k, n - k, deriv = deriv, integrate = integrate) /
               choose(n - k, j - k)
    ret * sqrt(2 * (n - j) + 1)
}



### set-up Bernstein polynom basis functions
### http://en.wikipedia.org/wiki/Bernstein_polynomial
### http://dx.doi.org/10.1080/02664761003692423
### arXiv preprint arXiv:1404.2293
Bernstein_basis <- function(order = 2, support = c(0, 1), normal = FALSE,
                            ui = c("none", "increasing", "decreasing", "cyclic", "zerointegral"), 
                            ci = 0, varname = NULL) {

    zeroint <- FALSE
    ui <- match.arg(ui, several.ok = TRUE)
    if (length(ui) > 2) ui <- "none"
    if (length(ui) > 1) 
        stopifnot((length(ui) == 2) && "zerointegral" %in% ui)
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

        if (!normal || ui != "none")
            X <- do.call("cbind", lapply(0:order, function(j) 
                         .Bx(x, j, order, deriv = deriv, integrate = integrate)))
        if (normal) {
            oX <- do.call("cbind", lapply(0:order, function(j) 
                          .oBx(x, j, order, deriv = deriv, integrate = integrate)))
            if (ui != "none" & (deriv == 0 & !integrate)) {
                ### transform constraint on raw Bernstein bases into
                ### constraint on orthogonal bases
                XtX <- crossprod(X)
                ### did not work well oui <- try(constr$ui %*% chol2inv(chol(XtX)) %*% crossprod(X, oX))
                oui <- try(constr$ui %*% solve(XtX) %*% crossprod(X, oX))
                oci <- constr$ci
                if (inherits(oui, "try-error")) { ### use f'(x) >= 0
                    #ngrid <- max(c(2 * order, 20))
                    #xgrid <- 0:ngrid / ngrid
                    #oui <- do.call("cbind", lapply(0:order, function(j)
                    #    .oBx(xgrid, j, order, deriv = 1, integrate = FALSE)))
                    #if (ui == "decreasing") oui <- -oui
                    #oci <- rep(constr$ci, length = nrow(oui))
                    warning("cannot obtain constraints for orthogonal Bernstein basis, 
                            returning unconstraint basis")
                    oX <- X
                    oui <- constr$ui
                }
                constr$ui <- oui
                constr$ci <- oci
            }
            X <- oX
        }
        if (deriv > 0)
            X <- X * (1 / diff(support)^deriv)
        colnames(X) <- paste("Bs", 1:ncol(X), "(", varname, ")", sep = "")
        if (zeroint) { ### normal???
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
    attr(basis, "intercept") <- zeroint

    class(basis) <- c("Bernstein_basis", "basis", class(basis))
    return(basis)
}

### evaluate model.matrix of Bernstein polynom
model.matrix.Bernstein_basis <- function(object, data,
                                         deriv = 0L, integrate = FALSE, ...)
    model.matrix.basis(object = object, data = data, 
                       deriv = deriv, integrate = integrate, ...)
