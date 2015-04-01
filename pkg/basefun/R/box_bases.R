
### box product of multiple _named_ basi(e)s objects
b <- function(..., sumconstr = FALSE) {

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = c("basis", "bases"))))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bases))

    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))

    attr(bases, "sumconstr") <- sumconstr

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
    if (attr(object, "sumconstr")) {
        s1 <- generate(object[[1]], 2)
        X1 <- model.matrix(object[[1]], data = expand.grid(s1))
        s2 <- generate(object[[2]], 2) ### min/max for numerics; all levels for factors
        X2 <- model.matrix(object[[2]], data = expand.grid(s2))
        ui <- attr(X1, "constraint")$ui
        ci <- attr(X1, "constraint")$ci
        ui <- kronecker(X2, ui)
        ci <- rep(ci, nrow(X2))
        constr <- list(ui = ui, ci = ci)
    } else {
        constr <- do.call(".box_ui_ci", lapply(ret, function(r)
                          attr(r, "constraint")))
    }
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
    vn <- lapply(object, varnames)
#    stopifnot(all(sapply(vn, length) == 1))
    stopifnot(all(!duplicated(unlist(vn))))
    X <- model.matrix(object = object, data = newdata, 
                      model.matrix = FALSE, ...)
    X <- lapply(X, function(x) as(x, "matrix"))
    X$beta <- array(coef, sapply(X, NCOL))
    ret <- do.call(".cXb", X)
#    dimnames(ret) <- names(object)
    ret
}
