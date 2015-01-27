
factor_basis <- function(support, contrast = "contr.treatment", remove_intercept = TRUE,
                         constraint = c("none"), varname = NULL) {

    constraint <- match.arg(constraint)

    basis <- function(x) {
        stopifnot(all(x %in% support))
        stopifnot(isTRUE(all.equal(class(x), class(support))))
        X <- model.matrix(~ x, data = data.frame(x = x), contrasts.arg = list(x = contrast))
        if (remove_intercept)
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        attr(X, "constraint") <- switch(constraint,
            "none" = list(ui = Diagonal(ncol(X)), 
                      ci = rep(-Inf, ncol(X))))
        return(X)
    }

    attr(basis, "length") <- nlevels(support) - remove_intercept
    attr(basis, "support") <- support
    attr(basis, "dimension") <- NULL
    attr(basis, "varnames") <- varname


    class(basis) <- c("basis", "factor_basis", class(basis))
    return(basis)
}

