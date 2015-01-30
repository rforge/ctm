
as.bases <- function(object, ...)
    UseMethod("as.bases")

as.bases.formula <- function(object, remove_intercept = FALSE, 
                             vars = NULL, varnames = all.vars(object),
                             constraint = NULL, ...) {

    ret <- function(data, deriv = 0, ...) {
        data <- data[varnames]
        data <- as.data.frame(data)
        if (deriv == 1 & length(varnames) == 1) {
            data[[varnames]] <- 1
        } else {
            deriv <- 0
        }
        mf <- model.frame(object, data, ...)
        X <- model.matrix(object, data = mf, ...)
        if (remove_intercept) {
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        } else if (deriv == 1) {
            X[, colnames(X) == "(Intercept)"] <- 0
        }
        if (is.null(constraint))
            constraint <- list(ui = Diagonal(ncol(X)),
                               ci = rep(-Inf, ncol(X)))
        attr(X, "constraint") <- constraint
        return(X)
    }
    nc <- NA
    s <- NA
    if (!is.null(vars)) {
        Xtmp <- ret(vars)
        nc <- ncol(Xtmp)
        s <- lapply(vars, function(v) {
            if (is.factor(v)) return(unique(v))
            return(NULL)
        })
    }
    attr(ret, "length") <- nc
    attr(ret, "varnames") <- varnames
    attr(ret, "dimension") <- nc
    attr(ret, "support") <- s
    class(ret) <- c("bases", "basis", class(ret))
    tmp <- ret
    class(tmp) <- class(tmp)[-1]
    bases <- list(self = tmp)
    ret
}
