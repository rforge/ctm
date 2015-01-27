

### set-up Bernstein polynom basis functions
### http://en.wikipedia.org/wiki/Bernstein_polynomial
### http://dx.doi.org/10.1080/02664761003692423
Bernstein_basis <- function(order = 2, support = c(0, 1),
                            constraint = c("none", "increasing", "decreasing"),
                            varname = NULL) {

    constraint <- match.arg(constraint)
    constr <- switch(constraint,
        "none" = list(ui = Diagonal(order + 1), 
                      ci = rep(-Inf, order + 1)),
        "increasing" = list(ui = diff(Diagonal(order + 1), differences = 1), 
                            ci = rep(0, order)),
        "decreasing" = list(ui = diff(Diagonal(order + 1), differences = 1) * -1,
                            ci = rep(0, order)))

    basis <- function(x, deriv = 0L, integrate = FALSE) {
        stopifnot(order > deriv)
        x <- (x - support[1]) / diff(support)
        fun <- ifelse(isTRUE(integrate), pbeta, dbeta)
        X <- do.call("cbind", lapply(0:(order - deriv), 
            function(m) fun(x, shape1 = m + 1, 
                            shape2 = (order - deriv) - m + 1) / ((order - deriv) + 1)))
        if (deriv > 0L) {
            fact <- prod(order:(order - deriv + 1)) * (1 / diff(support)^deriv)
            X <- X %*% diff(diag(order + 1), differences = deriv) * fact
        }
        attr(X, "constraint") <- constr
        return(X)
    }

    attr(basis, "length") <- order + 1L
    attr(basis, "dimension") <- NULL
    attr(basis, "support") <- support
    attr(basis, "varnames") <- varname

    class(basis) <- c("Bernstein_basis", "basis", class(basis))
    return(basis)
}

### evaluate model.matrix of Bernstein polynom
model.matrix.Bernstein_basis <- function(object, data, 
                                         deriv = 0L, integrate = FALSE, ...)
    model.matrix.basis(object = object, data = data, 
                       deriv = deriv, integrate = integrate, ...)

