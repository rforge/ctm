
### box product of multiple _named_ basi(e)s objects
b <- function(...) {

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = "basis")))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bnames))

    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))

    class(bases) <- c("box_bases", "bases")
    bases
}

model.matrix.box_bases <- function(object, data, model.matrix = TRUE,
    deriv = NULL, integrate = NULL, ...) {

    varnames <- varnames(object)
    bnames <- names(object)

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
            if (names(deriv) %in% varnames(object[[b]]))
                thisargs$deriv <- deriv
        }
        if (!is.null(integrate)) {
            if (names(integrate) %in% varnames(object[[b]]))
                thisargs$integrate <- integrate
        }
        thisargs$object <- object[[b]]
        thisargs$data <- data
        X <- do.call("model.matrix", thisargs)
        attr(X, "Assign") <- rbind(attr(X, "Assign"), b)
        X
    })
    if (!model.matrix) return(ret)
    constr <- do.call(".box_ui_ci", lapply(ret, function(r)
                      attr(r, "constraint")))
    if (length(object) > 1) {
        a <- do.call(".box_char", lapply(ret, function(r) attr(r, "Assign")))
        ret <- do.call(".box", ret)
        attr(ret, "Assign") <- a
    } else {
        ret <- ret[[1]]
    }
    attr(ret, "constraint") <- constr
    return(ret )
}

nparm.box_bases <- function(object, data)
    prod(sapply(object, nparm, data = data))

predict.box_bases <- function(object, newdata, coef, ...) {
    if (is.data.frame(newdata))    
        return(predict.basis(object = object, newdata = newdata,
                             coef = coef, ...))
    X <- model.matrix(object = object, data = newdata, 
                      model.matrix = FALSE, ...)
    X <- lapply(X, function(x) as(x, "matrix"))
    X$beta <- array(coef, sapply(X, NCOL))
    do.call(".cXb", X)
}
