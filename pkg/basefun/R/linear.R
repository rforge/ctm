
### linear basis
linear_basis <- function(support = c(0, 1), remove_intercept = FALSE,
                         constraint = c("none", "increasing", "decreasing"),
                         varname = NULL) {

    constraint <- match.arg(constraint)
    constr <- switch(constraint,
        "none" = list(ui = Diagonal(2), 
                      ci = rep(-Inf, -Inf)),
        "increasing" = list(ui = Diagonal(2),
                            ci = c(-Inf, 0)),
        "decreasing" = list(ui = Diagonal(2) * -1,
                            ci = c(-Inf, 0)))

    basis <- function(x) {
        x <- (x - support[1]) / diff(support)
        X <- cbind("(Intercept)" = 1, x = x)
        if (remove_intercept) {
            X <- X[,-1, drop = FALSE]
            constr$ui <- constr$ui[,-1, drop = FALSE]
            constr$ci <- constr$ci[-1]
        }
        attr(X, "constraint") <- constr
        return(X)
    }

    attr(basis, "length") <- 1 + remove_intercept
    attr(basis, "dimension") <- NULL
    attr(basis, "support") <- support
    attr(basis, "varnames") <- varname

    class(basis) <- c("linear_basis", "basis", class(basis))
    return(basis)
}
