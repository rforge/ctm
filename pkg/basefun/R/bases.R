
### evaluate model.matrix of multiple basis objects
model.matrix.bases <- function(object, data, ...)
    object(data, ...)

### evaluate model.matrix of single basis object
model.matrix.basis <- function(object, data, ...) {
    vars <- varnames(object)
    if (is.null(vars)) vars <- 1
    object(data[[vars]], ...)
}

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
    X$beta <- array(coef, dim(object))
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

length.basis <- function(x)
    attr(x, "length")

dim.basis <- function(x)
    attr(x, "dimension")

### box product of multiple _named_ basi(e)s objects
b <- function(...) {

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = "basis")))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bnames))

    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))
    support <- .rec2flat(lapply(bases, support))
    names(support) <- varnames ### is this safe???
    length <- prod(sapply(bases, length))
    dim <- sapply(bases, length)

    ret <- function(data, model.matrix = TRUE, ...) {
        args <- list(...)
        stopifnot(all(names(args) %in% bnames))
        ret <- lapply(bnames, function(b) {
            thisargs <- args[[b]]
            thisargs$object <- bases[[b]]
            thisargs$data <- data
            do.call("model.matrix", thisargs)
        })
        if (!model.matrix) return(ret)
        constr <- do.call(".box_ui_ci", lapply(ret, function(r)
                          attr(r, "constraint")))
        if (length(bases) > 1) {
            ret <- do.call(".box", ret)
        } else {
            ret <- ret[[1]]
        }
        attr(ret, "constraint") <- constr
        return(ret )
    }
    attr(ret, "length") <- length
    attr(ret, "varnames") <- varnames
    attr(ret, "dimension") <- dim
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
    support <- .rec2flat(lapply(bases, support))
    names(support) <- varnames ### is this safe???
    length <- prod(sapply(bases, length))
    dim <- sapply(bases, length)

    ret <- function(data, model.matrix = TRUE, ...) {
        args <- list(...)
        stopifnot(all(names(args) %in% bnames))
        ret <- lapply(bnames, function(b) {
            thisargs <- args[[b]]
            thisargs$object <- bases[[b]]
            thisargs$data <- data
            do.call("model.matrix", thisargs)
        })
        if (!model.matrix) return(ret)
        ui <- do.call("bdiag", lapply(ret, function(r)
                      attr(r, "constraint")$ui))
        ci <- do.call("c", lapply(ret, function(r)
                      attr(r, "constraint")$ci))
        if (length(bases) > 1) {
            ret <- do.call("cbind", ret)
        } else {
            ret <- ret[[1]]
        }
        attr(ret, "constraint") <- list(ui = ui, ci = ci)
        return(ret )
    }
    attr(ret, "length") <- length
    attr(ret, "varnames") <- varnames
    attr(ret, "dimension") <- dim
    attr(ret, "support") <- support
    class(ret) <- c("cbases", "bases", "basis", class(ret))
    ret
}

generate <- function(object, n)
    UseMethod("generate")

.generate <- function(s, n) {
    if ((length(s) == 2) && (storage.mode(s) == "double")) {
        x <- seq(from = s[1], to = s[2], length.out = n)
    } else if (is.factor(s)) {
        x <- s
    } else {
        x <- NA
    }
    return(x)
}

generate.basis <- function(object, n) {
    .generate(support(object), n = n)
}
 
generate.bases <- function(object, n) {
    s <- support(object)
    ret <- lapply(s, .generate, n = n)
    ret[!duplicated(names(ret))]
}

