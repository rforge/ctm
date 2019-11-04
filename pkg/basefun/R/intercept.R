
intercept_basis <- function(ui = c("none", "increasing", "decreasing"), negative = FALSE) {

    ui <- match.arg(ui)
    ci <- switch(ui, "none" = -Inf, 0)
    ui <- switch(ui, "none" = matrix(1),
                     "increasing" = matrix(1),
                     "decreasing" = matrix(-1))

    basis <- function(data, deriv = 0L) {

        X <- matrix(1, nrow = NROW(data))
        if (negative) X <- -X
        if (deriv > 0L) {
            if (!"(_Intercept_)" %in% names(deriv)) X[] <- 0
        }
        if (deriv < 0L) X[] <- 0
        colnames(X) <- "(_Intercept_)"
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- "(_Intercept_)"
        X
    }

    attr(basis, "variables") <- NA
    attr(basis, "intercept") <- TRUE

    class(basis) <- c("intercept_basis", "basis", class(basis))
    return(basis)
}

model.matrix.intercept_basis <- function(object, data, deriv = 0L, 
                                         dim = NULL, ...) {
    ### this is an intercept, so coef * 1 is everything we need
    if (!is.null(dim)) data <- data.frame(1)
    object(data, deriv = .deriv("(_Intercept_)", deriv))
}
