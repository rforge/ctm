
as.basis.formula <- function(object, remove_intercept = FALSE, 
                             data = NULL, ci = NULL, ui = NULL, ...) {

    varnames <- all.vars(object)

    mf <- NULL
    if (!is.null(data) && NROW(data) > 0) {
        ### see http://developer.r-project.org/model-fitting-functions.txt
        mf <- model.frame(object, data = data, drop.unused.levels = TRUE)
        mt <- attr(mf, "terms")
        X <- model.matrix(mt, data = mf, ...)
        contr <- attr(X, "contrasts")
        xlevels <- .getXlevels(mt, mf)
    }
    if (!is.null(data)) {
        s <- lapply(data[, varnames, drop = FALSE], function(v) {
            if (is.factor(v)) return(unique(v))
            if (length(v) > 0)
                return(range(v, na.rm = TRUE))
            return(NULL)
        })
    } else {
        s <- NULL
    }

    ret <- function(data) {
        data <- data[varnames]
        stopifnot(is.data.frame(data)) 
        if (!is.null(mf)) {
            mf <- model.frame(mt, data = data, xlev = xlevels)
            if(!is.null(cl <- attr(mt, "dataClasses"))) .checkMFClasses(cl, mf)
            X <- model.matrix(mt, mf, contrasts = contr)
        } else {
            mf <- model.frame(object, data)
            X <- model.matrix(attr(mf, "terms"), data = mf, ...)
        }
        if (remove_intercept) {
            a <- attr(X, "assign")
            X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
            attr(X, "assign") <- a[-1]
        }
        if (is.null(ui)) ui <- Diagonal(ncol(X))
        if (is.null(ci)) ci <- rep(-Inf, ncol(X))
        attr(X, "constraint") <- list(ui = ui, ci = ci)
        attr(X, "Assign") <- c("(Intercept)", varnames)[attr(X, "assign") + 1]
        return(X)
    }

    attr(ret, "varnames") <- varnames
    attr(ret, "support") <- s
    ### note: ~ a - 1 also contains intercept!
    attr(ret, "intercept") <- !remove_intercept 
    class(ret) <- c("formula_basis", "basis", class(ret))
    ret
}

model.matrix.formula_basis <- function(object, data, dim = NULL, ...) {

    if (!is.null(dim)) {
        nd <- names(dim)
        data <- do.call("expand.grid", data[nd[nd %in% varnames(object)]])
    }

    object(data)
}
