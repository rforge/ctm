
### evaluate model.matrix of multiple basis objects
model.matrix.basis <- function(object, data, ...)
    object(data, ...)

### compute predicted values of (multiple) basis objects
predict.basis <- function(object, newdata, coef, ...) {
    X <- model.matrix(object = object, data = newdata, ...)
    X %*% coef
}

predict.bbases <- function(object, newdata, coef, ...) {
    if (is.data.frame(newdata)) 
        return(predict.basis(object = object, newdata = newdata, 
                             coef = coef, ...))
    X <- model.matrix(object = object, data = newdata, model.matrix = FALSE, ...)
    X <- lapply(X, function(x) as(x, "matrix"))
    X$beta <- array(coef, sapply(X, NCOL))
    do.call("cXb", X)
}

predict.cbases <- function(object, newdata, coef, ...) {
    stopifnot(is.data.frame(newdata)) 
    predict.basis(object = object, newdata = newdata, 
                  coef = coef, ...)
}

varnames <- function(x)
    UseMethod("varnames")

varnames.default <- function(x)
    attr(x, "varnames")

support <- function(x)
    UseMethod("support")

support.default <- function(x)
    attr(x, "support")

### box product of multiple _named_ basi(e)s objects
b <- function(...) {

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = "basis")))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bnames))

    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))
    support <- lapply(bases, support)

    ret <- function(data, model.matrix = TRUE, deriv = NULL, integrate = NULL) {
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
                if (names(deriv) %in% varnames(bases[[b]]))
                    thisargs$deriv <- deriv
            }
            if (!is.null(integrate)) {
                if (names(integrate) %in% varnames(bases[[b]]))
                    thisargs$integrate <- integrate
            }
            thisargs$object <- bases[[b]]
            thisargs$data <- data
            X <- do.call("model.matrix", thisargs)
            attr(X, "Assign") <- rbind(attr(X, "Assign"), b)
            X
        })
        if (!model.matrix) return(ret)
        constr <- do.call(".box_ui_ci", lapply(ret, function(r)
                          attr(r, "constraint")))
        if (length(bases) > 1) {
            a <- do.call(".box_char", lapply(ret, function(r) attr(r, "Assign")))
            ret <- do.call(".box", ret)
            attr(ret, "Assign") <- a
        } else {
            ret <- ret[[1]]
        }
        attr(ret, "constraint") <- constr
        return(ret )
    }
    attr(ret, "varnames") <- varnames
    attr(ret, "support") <- support
    class(ret) <- c("bbases", "bases", "basis", class(ret))
    ret
}

### cbind of multiple _named_ basi(e)s objects
c.basis <- function(..., recursive = FALSE) {

    stopifnot(!recursive)

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = "basis")))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bnames))

    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))
    support <- lapply(bases, support)

    ret <- function(data, model.matrix = TRUE, deriv = NULL, integrate = NULL) {
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
                if (names(deriv) %in% varnames(bases[[b]])) {
                    thisargs$deriv <- deriv
                } else {
                    X <- model.matrix(bases[[b]], data)
                    X[] <- 0
                    return(X)
                }
            }
            if (!is.null(integrate)) {
                if (names(integrate) %in% varnames(bases[[b]]))
                    thisargs$integrate <- integrate
            }
            thisargs$object <- bases[[b]]
            thisargs$data <- data
            X <- do.call("model.matrix", thisargs)
            X <- do.call("model.matrix", thisargs)
            attr(X, "Assign") <- rbind(attr(X, "Assign"), b)
            X
        })
        if (!model.matrix) return(ret)
        ui <- do.call("bdiag", lapply(ret, function(r)
                      attr(r, "constraint")$ui))
        ci <- do.call("c", lapply(ret, function(r)
                      attr(r, "constraint")$ci))
        if (length(bases) > 1) {
            a <- lapply(ret, function(r) matrix(attr(r, "Assign"), ncol = ncol(r)))
            mr <- max(sapply(a, NROW))
            for (i in 1:length(a)) a[[i]] <- a[[i]][rep(1:nrow(a[[i]]), length = mr),,drop = FALSE]
            ret <- do.call("cbind", ret)
            attr(ret, "Assign") <- do.call("cbind", a)
        } else {
            ret <- ret[[1]]
        }
        attr(ret, "constraint") <- list(ui = ui, ci = ci)
        return(ret )
    }
    attr(ret, "varnames") <- varnames
    attr(ret, "support") <- support
    class(ret) <- c("cbases", "bases", "basis", class(ret))
    ret
}

generate <- function(object, n)
    UseMethod("generate")


generate.basis <- function(object, n) {
    ret <- list()
    .generate <- function(s, n) {
        if (!is.atomic(s)) {
            tmp <- lapply(s, .generate, n = n)
            if (all(sapply(s, is.atomic)))
                ret[names(tmp)] <<- tmp
        }
        if ((length(s) == 2) && (storage.mode(s) == "double")) {
            x <- seq(from = s[1], to = s[2], length.out = n)
        } else if (is.factor(s)) {
            x <- s
        } else {
            x <- NULL
        }
        return(x)
    }
    ret2 <- .generate(support(object), n = n)
    if (length(ret) == 0) return(ret2)
    ret
}

