
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

nparm.cbind_bases <- function(object, data)
    sum(sapply(object, nparm, data = data))

model.matrix.cbind_bases <- function(object, data, model.matrix = TRUE, 
    deriv = NULL, integrate = NULL, ...) {

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

predict.cbind_bases <- function(object, newdata, coef, deriv = NULL, 
                                integrate = NULL, ...) {

    if (is.data.frame(newdata))
        return(predict.basis(object = object, newdata = newdata,
                             coef = coef, ...))

#    np <- sapply(object, nparm, data = newdata)
    if (!all(names(object) %in% names(newdata))) {
        newdata <- lapply(names(object), function(b) {
                          newdata[varnames(object[[b]])]
        })
        names(newdata) <- names(object)
    }
    np <- sapply(names(object), function(b) nparm(object[[b]], newdata[[b]]))

    ret <- lapply(1:length(object), function(b) {
        start <- ifelse(b == 1, 1, sum(np[1:(b - 1)]) + 1)
        cf <-  coef[start:sum(np[1:b])]
        thisargs <- list()
        if (!is.null(deriv)) { 
            if (names(deriv) %in% varnames(object[[b]])) {
                thisargs$deriv <- deriv
            } else {
                return(rep(0, nrow(model.matrix(object[[b]], data = newdata))))
	    }
        }
        if (!is.null(integrate)) {
            if (names(integrate) %in% varnames(object[[b]]))
                thisargs$integrate <- integrate
        }
        thisargs$object <- object[[b]]
        thisargs$newdata <- newdata[[names(object)[b]]]
        thisargs$coef <- cf
        return(do.call("predict", thisargs))
    })
    vn <- sapply(object, varnames)
    if (length(unique(varnames(object))) == 1) return(Reduce("+", ret))
#    stopifnot(all(sapply(vn, length) == 1))
    stopifnot(all(!duplicated(unlist(vn))))
    return(array(rowSums(expand.grid(ret)), dim = sapply(ret, NROW))) ### CHECK, outer?
}
