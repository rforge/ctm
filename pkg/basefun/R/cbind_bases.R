
### cbind of multiple _named_ basi(e)s objects
c.basis <- function(..., recursive = FALSE) {

    stopifnot(!recursive)

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = c("basis", "bases"))))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bases))
    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))
    inter <- sapply(bases, intercept)
    if (sum(inter) > 1) 
        warning("more than one basis contains an intercept term")
 
    class(bases) <- c("cbind_bases", "bases")
    bases
}

c.bases <- c.basis

model.matrix.cbind_bases <- function(object, data, model.matrix = TRUE, 
    dim = NULL, deriv = NULL, integrate = NULL, ...) {

    if (model.matrix) stopifnot(is.data.frame(data))

    bnames <- names(object)
    varnames <- varnames(object)

    if (!is.null(deriv)) {
        stopifnot(length(deriv) == 1)
        if (!names(deriv) %in% varnames) deriv <- NULL
    }
    if (!is.null(integrate)) {
        stopifnot(length(integrate) == 1)
        if (!names(integrate) %in% varnames) integrate <- NULL
    }
    ret <- lapply(bnames, function(b) {
        thisargs <- list()
        if (!is.null(deriv)) {
            if (names(deriv) %in% varnames(object[[b]])) {
                thisargs$deriv <- deriv
            } else {
                X <- model.matrix(object[[b]], data)
                X[] <- 0
                return(X)
            }
        }
        if (!is.null(integrate)) {
            if (names(integrate) %in% varnames(object[[b]]))
                thisargs$integrate <- integrate
        }
        thisargs$object <- object[[b]]
        thisargs$data <- data
        if (!is.null(dim))
            thisargs$dim <- dim[names(dim) %in% varnames(object[[b]])]
        X <- do.call("model.matrix", thisargs)
        attr(X, "Assign") <- rbind(attr(X, "Assign"), b)
        X
    })
    if (!model.matrix) return(ret)
    ui <- do.call("bdiag", lapply(ret, function(r)
                  attr(r, "constraint")$ui))
    ci <- do.call("c", lapply(ret, function(r)
                  attr(r, "constraint")$ci))
    if (length(object) > 1) {
        a <- lapply(ret, function(r) matrix(attr(r, "Assign"), 
                                            ncol = ncol(r)))
        mr <- max(sapply(a, NROW))
        for (i in 1:length(a)) 
            a[[i]] <- a[[i]][rep(1:nrow(a[[i]]), length = mr),,
                             drop = FALSE]
        ret <- do.call("cbind", ret)
        attr(ret, "Assign") <- do.call("cbind", a)
    } else {
        ret <- ret[[1]]
    }
    attr(ret, "constraint") <- list(ui = ui, ci = ci)
    return(ret )
}

predict.cbind_bases <- function(object, newdata, coef, 
                                dim = !is.data.frame(newdata), 
                                terms = names(object), ...) {

    if (isTRUE(dim))
        dim <- sapply(newdata, NROW) 
    else if (is.logical(dim)) 
        dim <- FALSE

    if (length(object) == 1) 
        return(predict(object[[1]], newdata = newdata, 
                       coef = coef, dim = dim, ...))

    np <- nparm(object, data = newdata)

    ret <- vector(mode = "list", length = length(object))
    names(ret) <- names(object)
    for (b in 1:length(object)) {
        if (names(object)[b] %in% terms || names(object[[b]]) %in% terms) {
            start <- ifelse(b == 1, 1, sum(unlist(np[names(object)[1:(b - 1)]])) + 1)
            cf <- coef[start:sum(unlist(np[names(object)[1:b]]))]
            ret[[b]] <- predict(object[[b]], newdata = newdata, 
                                coef = cf, dim = dim, terms = terms, ...)
        } else {
            ret[[b]] <- array(0, dim = dim)
        }
    }
    return(Reduce("+", ret))
}
